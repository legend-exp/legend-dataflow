# ruff: noqa: T201, I002

"""Run-completion script executed by the ``autogen_output`` rule: checks log
files for warnings and errors, builds ``pygama`` FileDB databases and
collects the produced run information."""

import datetime
import json
import os
import subprocess
import time
from pathlib import Path

from legenddataflowscripts.workflow import as_ro, execenv_pyexe
from legenddataflowscripts.workflow.execenv import _execenv2str

from legenddataflow.methods import paths as pat
from legenddataflow.methods import patterns

try:
    from snakemake.script import snakemake
except ImportError:  # not running inside a snakemake job (e.g. unit tests)
    snakemake = None


def check_log_files(log_path, output_file, gen_output, warning_file=None):
    now = datetime.datetime.now(datetime.UTC).strftime("%d/%m/%y %H:%M")

    # collect ERROR/WARNING lines from all log files, consuming the logs;
    # uncaught tracebacks carry no "ERROR" token, so track traceback blocks
    # and report their final exception line as an error
    error_lines = []
    warning_lines = []
    for file in Path(log_path).rglob("*.log"):
        in_traceback = False
        with Path(file).open() as r:
            for line in r:
                if line.startswith("Traceback (most recent call last):"):
                    in_traceback = True
                elif in_traceback:
                    if line[:1] in (" ", "\t") or not line.strip():
                        continue
                    # first unindented line is the exception itself
                    error_lines.append(f"{Path(file).name} : {line.rstrip()}\n")
                    in_traceback = False
                elif "ERROR" in line:
                    error_lines.append(f"{Path(file).name} : {line.rstrip()}\n")
                elif "WARNING" in line:
                    warning_lines.append(f"{Path(file).name} : {line.rstrip()}\n")
        Path(file).unlink()

    def _write_summary(path, lines, label):
        Path(path).parent.mkdir(parents=True, exist_ok=True)
        with Path(path).open("w") as f:
            if lines:
                f.write(f"{gen_output} successfully generated at {now} with {label} \n")
                f.writelines(lines)
            else:
                f.write(
                    f"{gen_output} successfully generated at {now} with no {label} \n"
                )

    _write_summary(output_file, error_lines, "errors")
    if warning_file is not None:
        _write_summary(warning_file, warning_lines, "warnings")

    walk = list(os.walk(log_path))
    for path, _, _ in walk[::-1]:
        if len(list(Path(path).iterdir())) == 0:
            Path(path).rmdir()


def find_gen_runs(gen_tier_path):
    # first look for non-concat tiers
    paths = gen_tier_path.glob("*/*/*/*")
    # use the directories to build a datatype/period/run string
    runs = {"/".join(str(p).split("/")[-3:]) for p in paths}

    # then look for concat tiers (use filenames now)
    paths_concat = gen_tier_path.glob("*/*/*.lh5")
    # {experiment}-{period}-{run}-{datatype}-... -> datatype/period/run
    runs_concat = {
        "/".join([p.name.split("-")[3], *p.name.split("-")[1:3]]) for p in paths_concat
    }

    return runs | runs_concat


def build_file_dbs(setup, gen_tier_path, outdir, ignore_keys_file=None, threads=1):
    tic = time.time()

    gen_tier_path = Path(as_ro(setup, gen_tier_path))
    outdir = Path(outdir)

    # find generated directories
    runs = find_gen_runs(gen_tier_path)

    if not runs:
        print(f"WARNING: did not find any processed runs in {gen_tier_path}")

    processes = set()
    for spec in runs:
        speck = spec.split("/")
        outdir.mkdir(parents=True, exist_ok=True)
        # TODO: replace l200 with {experiment}
        outfile = outdir / f"l200-{speck[1]}-{speck[2]}-{speck[0]}-filedb.h5"
        logfile = (
            Path(pat.tmp_log_path(setup)) / "filedb" / outfile.with_suffix(".log").name
        )

        print(f"INFO: ......building {outfile}")
        pre_cmdline, cmdenv = execenv_pyexe(setup, "build-filedb", as_string=False)

        cmdline = [
            *pre_cmdline,
            "--scan-path",
            spec,
            "--output",
            str(outfile),
            "--config",
            str(outdir / "file_db_config.json"),
            "--log",
            str(logfile),
        ]
        if ignore_keys_file is not None:
            cmdline += ["--ignore-keys", str(ignore_keys_file)]

        if speck[0] == "phy":
            cmdline += ["--assume-nonsparse"]

        # TODO: forward stdout to log file
        processes.add(subprocess.Popen(cmdline, env=cmdenv))

        if len(processes) >= threads:
            os.wait()
            processes.difference_update([p for p in processes if p.poll() is not None])

    for p in processes:
        if p.poll() is None:
            p.wait()

    for p in processes:
        if p.returncode != 0:
            msg = f"at least one FileDB building thread failed: {_execenv2str(p.args, cmdenv)}"
            raise RuntimeError(msg)

    toc = time.time()
    dt = datetime.timedelta(seconds=toc - tic)
    print(f"INFO: ...took {dt}")


def fformat(setup, tier):
    abs_path = patterns.get_pattern_tier(setup, tier, check_in_cycle=False)
    return str(abs_path).replace(pat.get_tier_path(setup, tier), "")


def build_file_db_config(setup, filedb_path):
    file_db_config = {}
    prodenv_var = os.getenv("PRODENV")

    if prodenv_var is not None and prodenv_var in str(filedb_path):
        prodenv = as_ro(setup, prodenv_var)

        def tdirs(tier):
            return as_ro(setup, pat.get_tier_path(setup, tier)).replace(prodenv, "")

        file_db_config["data_dir"] = "$PRODENV"

    else:
        print("WARNING: $PRODENV not set, the FileDB will not be relocatable")

        def tdirs(tier):
            return as_ro(setup, pat.get_tier_path(setup, tier))

        file_db_config["data_dir"] = "/"

    file_db_config["tier_dirs"] = {k: tdirs(k) for k in setup["table_format"]}

    file_db_config |= {
        "file_format": {k: fformat(setup, k) for k in setup["table_format"]},
        "table_format": setup["table_format"],
    }
    return file_db_config


if snakemake is not None:
    print("INFO: dataflow ran successfully, now few final checks and scripts")
    setup = snakemake.params.setup

    if (snakemake.wildcards.tier not in ("daq", "daq_compress")) and setup.get(
        "build_file_dbs", True
    ):
        print(f"INFO: ...building FileDBs with {snakemake.threads} threads")

        file_db_config = build_file_db_config(setup, snakemake.params.filedb_path)

        Path(snakemake.params.filedb_path).mkdir(parents=True, exist_ok=True)
        with (Path(snakemake.params.filedb_path) / "file_db_config.json").open(
            "w"
        ) as f:
            json.dump(file_db_config, f, indent=2)

        build_file_dbs(
            setup,
            pat.tier_path(setup),
            snakemake.params.filedb_path,
            snakemake.params.ignore_keys_file,
            threads=snakemake.threads,
        )
        (Path(snakemake.params.filedb_path) / "file_db_config.json").unlink()

    if setup.get("check_log_files", True):
        print("INFO: ...checking log files")
        check_log_files(
            pat.tmp_log_path(setup),
            snakemake.output.summary_log,
            snakemake.output.gen_output,
            warning_file=snakemake.output.warning_log,
        )
    else:
        Path(snakemake.output.summary_log).parent.mkdir(parents=True, exist_ok=True)
        Path(snakemake.output.summary_log).touch()
        Path(snakemake.output.warning_log).touch()

    Path(snakemake.output.gen_output).touch()
