# ruff: noqa: F821, T201

import datetime
import glob
import json
import os
import time
from pathlib import Path

import util.utils as ut
from util import patterns
from util.FileKey import FileKey

print("INFO: dataflow ran successfully, now few final checks and scripts")


def as_ro(path):
    return ut.as_ro(snakemake.params.setup, path)


def check_log_files(log_path, output_file, gen_output, warning_file=None):
    now = datetime.datetime.now(datetime.UTC).strftime("%d/%m/%y %H:%M")
    os.makedirs(os.path.dirname(output_file), exist_ok=True)
    if warning_file is not None:
        os.makedirs(os.path.dirname(warning_file), exist_ok=True)
        with open(warning_file, "w") as w, open(output_file, "w") as f:
            n_errors = 0
            n_warnings = 0
            for file in Path(log_path).rglob("*.log"):
                with open(file) as r:
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
                                f.write(f"{os.path.basename(file)} : {line}\n")
                                n_errors += 1
                            elif "WARNING" in line:
                                w.write(f"{os.path.basename(file)} : {line}\n")
                                n_warnings += 1
                    else:
                        pass
                os.remove(file)
                text = None
            if n_errors == 0:
                f.write(f"{gen_output} successfully generated at {now} with no errors \n")
            if n_warnings == 0:
                w.write(f"{gen_output} successfully generated at {now} with no warnings \n")
    else:
        with open(output_file, "w") as f:
            n_errors = 0
            for file in Path(log_path).rglob("*.log"):
                with open(file) as r:
                    text = r.read()
                    if "ERROR" in text:
                        for line in text.splitlines():
                            if "ERROR" in line:
                                if n_errors == 0:
                                    f.write(
                                        f"{gen_output} successfully generated at {now} with errors \n"
                                    )
                                f.write(f"{os.path.basename(file)} : {line}\n")
                                n_errors += 1
                    else:
                        pass
                os.remove(file)
                text = None
            if n_errors == 0:
                f.write(f"{gen_output} successfully generated at {now} with no errors \n")
    walk = list(os.walk(log_path))
    for path, _, _ in walk[::-1]:
        if len(os.listdir(path)) == 0:
            os.rmdir(path)


def add_spaces(n):
    out_string = ""
    for _i in range(n):
        out_string += " "
    return out_string


def readable_json(dic, ncol=6, indent=4):
    def reformat_dict(dic, out_string="", indent_level=0, ncol=ncol, indent=indent):
        for key, item in dic.items():
            if isinstance(item, list):
                out_string += (
                    f'{add_spaces(indent_level+indent)}"{key}"'
                    + ": [\n"
                    + f"{add_spaces(2*indent+indent_level)}"
                )
                for i, _item in enumerate(item):
                    if i > 0 and (i) % ncol == 0:
                        out_string += f"\n{add_spaces(2*indent+indent_level)}"
                    out_string += f'"{_item}", '
                out_string = out_string[:-2]
                out_string += "\n" + f"{add_spaces(indent+indent_level)}" + "],\n"

            elif isinstance(item, dict):
                out_string += f'{add_spaces(indent+indent_level)}"{key}": ' + "{\n"
                out_string = reformat_dict(
                    item, out_string, indent_level=indent_level + indent, ncol=ncol
                )
                out_string += "\n" + f"{add_spaces(indent_level+indent)}" + "},\n"
        return out_string[:-2]

    out_string = "{\n"
    out_string = reformat_dict(dic, out_string=out_string, ncol=6)
    out_string += "\n}\n"

    return out_string


def get_keys(files):
    def get_run(Filekey):
        return f"{Filekey.experiment}-{Filekey.period}-{Filekey.run}-{Filekey.datatype}"

    key_dict = {}
    for file in files:
        key = FileKey.get_filekey_from_filename(os.path.basename(file))
        if get_run(key) in key_dict:
            key_dict[get_run(key)].append(file)
        else:
            key_dict[get_run(key)] = [file]
    return key_dict


def build_valid_keys(input_files, output_dir):
    infiles = glob.glob(as_ro(input_files))
    key_dict = get_keys(infiles)

    for key in list(key_dict):
        dtype = key.split("-")[-1]
        out_file = os.path.join(output_dir, f'{key.replace(f"-{dtype}", "")}-valid_{dtype}.json')
        Path(os.path.dirname(out_file)).mkdir(parents=True, exist_ok=True)
        if os.path.isfile(out_file):
            out_dict = Props.read_from([out_file] + key_dict[key])
        else:
            out_dict = Props.read_from(key_dict[key])
        out_string = readable_json(out_dict)
        with open(out_file, "w") as w:
            w.write(out_string)

    for input_file in infiles:
        if os.path.isfile(input_file):
            os.remove(input_file)


def find_gen_runs(gen_tier_path):
    # first look for non-concat tiers
    paths = gen_tier_path.glob("*/*/*/*")
    # use the directories to build a datatype/period/run string
    runs = {"/".join(p.name.split("/")[-3:]) for p in paths}

    # then look for concat tiers (use filenames now)
    paths_concat = gen_tier_path.glob("*/*/*.lh5")
    # use the directories to build a datatype/period/run string
    runs_concat = {"/".join([p.name.split("-")[3]] + p.name.split("-")[1:3]) for p in paths_concat}

    return runs | runs_concat


def build_file_dbs(gen_tier_path, outdir):
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
        logfile = Path(ut.tmp_log_path(snakemake.params.setup)) / outfile.with_suffix(".log").name
        print(f"INFO: ......building {outfile}")

        cmdline = ut.runcmd(snakemake.params.setup, aslist=True)
        prodenv = as_ro(os.getenv("PRODENV"))
        cmdline += [f"--env=PRODENV={prodenv}"]

        # TODO: forward stdout to log file
        processes.add(
            subprocess.Popen(
                [
                    *cmdline,
                    "python3",
                    "-B",
                    f"{snakemake.params.basedir}/scripts/build_fdb.py",
                    "--scan-path",
                    spec,
                    "--output",
                    str(outfile),
                    "--config",
                    str(outdir / "file_db_config.json"),
                    "--log",
                    str(logfile),
                    "--assume-nonsparse" if speck[0] == "phy" else "",
                ],
            )
        )

        if len(processes) >= snakemake.threads:
            os.wait()
            processes.difference_update([p for p in processes if p.poll() is not None])

    for p in processes:
        if p.poll() is None:
            p.wait()

    for p in processes:
        if p.returncode != 0:
            _cmdline = " ".join([f"{k}={v}" for k, v in cmdenv.items()]) + " " + " ".join(p.args)
            msg = f"at least one FileDB building thread failed: {_cmdline}"
            raise RuntimeError(msg)

    toc = time.time()
    dt = datetime.timedelta(seconds=toc - tic)
    print(f"INFO: ...took {dt}")


file_db_config = {}

if os.getenv("PRODENV") is not None and os.getenv("PRODENV") in snakemake.params.filedb_path:
    prodenv = as_ro(os.getenv("PRODENV"))

    def tdirs(tier):
        return as_ro(ut.get_tier_path(snakemake.params.setup, tier)).replace(prodenv, "")

    file_db_config["data_dir"] = "$PRODENV"

else:
    print("WARNING: $PRODENV not set, the FileDB will not be relocatable")

    def tdirs(tier):
        return as_ro(ut.get_tier_path(snakemake.params.setup, tier))

    file_db_config["data_dir"] = "/"


file_db_config["tier_dirs"] = {k: tdirs(k) for k in snakemake.params.setup["table_format"]}


def fformat(tier):
    return as_ro(
        patterns.get_pattern_tier(snakemake.params.setup, tier, check_in_cycle=False)
    ).replace(as_ro(ut.get_tier_path(snakemake.params.setup, tier)), "")


file_db_config |= {
    "file_format": {k: fformat(k) for k in snakemake.params.setup["table_format"]},
    "table_format": snakemake.params.setup["table_format"],
}

if snakemake.wildcards.tier != "daq":
    print(f"INFO: ...building FileDBs with {snakemake.threads} threads")

    os.makedirs(snakemake.params.filedb_path, exist_ok=True)

    with open(os.path.join(snakemake.params.filedb_path, "file_db_config.json"), "w") as f:
        json.dump(file_db_config, f, indent=2)

    build_file_dbs(ut.tier_path(snakemake.params.setup), snakemake.params.filedb_path)
    os.remove(os.path.join(snakemake.params.filedb_path, "file_db_config.json"))

    build_valid_keys(
        os.path.join(ut.tmp_par_path(snakemake.params.setup), "*_db.json"),
        snakemake.params.valid_keys_path,
    )

print("INFO: ...checking log files")

check_log_files(
    ut.tmp_log_path(snakemake.params.setup),
    snakemake.output.summary_log,
    snakemake.output.gen_output,
    warning_file=snakemake.output.warning_log,
)

Path(snakemake.output.gen_output).touch()
