# ruff: noqa: F821, T201
from __future__ import annotations

import datetime
import json
import os
import subprocess
import time
from pathlib import Path

from snakemake.script import snakemake

from legenddataflow import paths as pat
from legenddataflow import patterns
from legenddataflow import utils as ut
from legenddataflow.execenv import _execenv2str, execenv_pyexe

print("INFO: dataflow ran successfully, now few final checks and scripts")


def as_ro(path):
    return ut.as_ro(snakemake.params.setup, path)


def check_log_files(log_path, output_file, gen_output, warning_file=None):
    now = datetime.datetime.now(datetime.UTC).strftime("%d/%m/%y %H:%M")
    Path(output_file).parent.mkdir(parents=True, exist_ok=True)
    if warning_file is not None:
        Path(warning_file).parent.mkdir(parents=True, exist_ok=True)
        with Path(warning_file).open("w") as w, Path(output_file).open("w") as f:
            n_errors = 0
            n_warnings = 0
            for file in Path(log_path).rglob("*.log"):
                with Path(file).open() as r:
                    text = r.read()
                    if "ERROR" in text or "WARNING" in text:
                        for line in text.splitlines():
                            if "ERROR" in line:
                                if n_errors == 0:
                                    f.write(
                                        f"{gen_output} successfully generated at {now} with errors \n"
                                    )
                                if n_warnings == 0:
                                    w.write(
                                        f"{gen_output} successfully generated at {now} with warnings \n"
                                    )
                                f.write(f"{Path(file).name} : {line}\n")
                                n_errors += 1
                            elif "WARNING" in line:
                                w.write(f"{Path(file).name} : {line}\n")
                                n_warnings += 1
                    else:
                        pass
                Path(file).unlink()
                text = None
            if n_errors == 0:
                f.write(
                    f"{gen_output} successfully generated at {now} with no errors \n"
                )
            if n_warnings == 0:
                w.write(
                    f"{gen_output} successfully generated at {now} with no warnings \n"
                )
    else:
        with Path(output_file).open("w") as f:
            n_errors = 0
            for file in Path(log_path).rglob("*.log"):
                with Path(file).open() as r:
                    text = r.read()
                    if "ERROR" in text:
                        for line in text.splitlines():
                            if "ERROR" in line:
                                if n_errors == 0:
                                    f.write(
                                        f"{gen_output} successfully generated at {now} with errors \n"
                                    )
                                f.write(f"{Path(file).name} : {line}\n")
                                n_errors += 1
                    else:
                        pass
                Path(file).unlink()
                text = None
            if n_errors == 0:
                f.write(
                    f"{gen_output} successfully generated at {now} with no errors \n"
                )
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
    # use the directories to build a datatype/period/run string
    runs_concat = {
        "/".join([str(p).split("-")[3]] + str(p).split("-")[1:3]) for p in paths_concat
    }

    return runs | runs_concat


def build_file_dbs(gen_tier_path, outdir, ignore_keys_file=None):
    tic = time.time()

    gen_tier_path = Path(as_ro(gen_tier_path))
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
            Path(pat.tmp_log_path(snakemake.params.setup))
            / "filedb"
            / outfile.with_suffix(".log").name
        )

        print(f"INFO: ......building {outfile}")
        pre_cmdline, cmdenv = execenv_pyexe(
            snakemake.params.setup, "build-filedb", as_string=False
        )

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

        if len(processes) >= snakemake.threads:
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


def fformat(tier):
    abs_path = patterns.get_pattern_tier(
        snakemake.params.setup, tier, check_in_cycle=False
    )
    return str(abs_path).replace(pat.get_tier_path(snakemake.params.setup, tier), "")


if snakemake.params.setup.get("build_file_dbs", True):
    file_db_config = {}

    if os.getenv("PRODENV") is not None and os.getenv("PRODENV") in str(
        snakemake.params.filedb_path
    ):
        prodenv = as_ro(os.getenv("PRODENV"))

        def tdirs(tier):
            return as_ro(pat.get_tier_path(snakemake.params.setup, tier)).replace(
                prodenv, ""
            )

        file_db_config["data_dir"] = "$PRODENV"

    else:
        print("WARNING: $PRODENV not set, the FileDB will not be relocatable")

        def tdirs(tier):
            return as_ro(pat.get_tier_path(snakemake.params.setup, tier))

        file_db_config["data_dir"] = "/"

    file_db_config["tier_dirs"] = {
        k: tdirs(k) for k in snakemake.params.setup["table_format"]
    }

    file_db_config |= {
        "file_format": {k: fformat(k) for k in snakemake.params.setup["table_format"]},
        "table_format": snakemake.params.setup["table_format"],
    }

if (snakemake.wildcards.tier != "daq") and snakemake.params.setup.get(
    "build_file_dbs", True
):
    print(f"INFO: ...building FileDBs with {snakemake.threads} threads")

    Path(snakemake.params.filedb_path).mkdir(parents=True, exist_ok=True)

    with (Path(snakemake.params.filedb_path) / "file_db_config.json").open("w") as f:
        json.dump(file_db_config, f, indent=2)

    build_file_dbs(
        pat.tier_path(snakemake.params.setup),
        snakemake.params.filedb_path,
        snakemake.params.ignore_keys_file,
    )
    (Path(snakemake.params.filedb_path) / "file_db_config.json").unlink()

if snakemake.params.setup.get("check_log_files", True):
    print("INFO: ...checking log files")
    check_log_files(
        pat.tmp_log_path(snakemake.params.setup),
        snakemake.output.summary_log,
        snakemake.output.gen_output,
        warning_file=snakemake.output.warning_log,
    )

Path(snakemake.output.gen_output).touch()
