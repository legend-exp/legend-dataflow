# ruff: noqa: F821, T201

import glob
import json
import os
import pathlib
from datetime import datetime
from pathlib import Path

import util.patterns as pat
import util.utils as ut
from util.CalibCatalog import Props
from util.FileKey import FileKey


def check_log_files(log_path, output_file, gen_output, warning_file=None):
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
                                        f"{gen_output} successfully generated at {datetime.utcnow().strftime('%d/%m/%y %H:%M')} with errors \n"
                                    )
                                if n_warnings == 0:
                                    w.write(
                                        f"{gen_output} successfully generated at {datetime.utcnow().strftime('%d/%m/%y %H:%M')} with warnings \n"
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
                f.write(
                    f"{gen_output} successfully generated at {datetime.utcnow().strftime('%d/%m/%y %H:%M')} with no errors \n"
                )
            if n_warnings == 0:
                w.write(
                    f"{gen_output} successfully generated at {datetime.utcnow().strftime('%d/%m/%y %H:%M')} with no warnings \n"
                )
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
                                        f"{gen_output} successfully generated at {datetime.utcnow().strftime('%d/%m/%y %H:%M')} with errors \n"
                                    )
                                f.write(f"{os.path.basename(file)} : {line}\n")
                                n_errors += 1
                    else:
                        pass
                os.remove(file)
                text = None
            if n_errors == 0:
                f.write(
                    f"{gen_output} successfully generated at {datetime.utcnow().strftime('%d/%m/%y %H:%M')} with no errors \n"
                )
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


def get_keys(input_files):
    def get_run(Filekey):
        return f"{Filekey.experiment}-{Filekey.period}-{Filekey.run}-{Filekey.datatype}"

    files = glob.glob(input_files)
    key_dict = {}
    for file in files:
        key = FileKey.get_filekey_from_filename(os.path.basename(file))
        if get_run(key) in key_dict:
            key_dict[get_run(key)].append(file)
        else:
            key_dict[get_run(key)] = [file]
    return key_dict


def build_valid_keys(input_files, output_dir):
    key_dict = get_keys(input_files)

    for key in list(key_dict):
        dtype = key.split("-")[-1]
        out_file = os.path.join(output_dir, f'{key.replace(f"-{dtype}", "")}-valid_{dtype}.json')
        pathlib.Path(os.path.dirname(out_file)).mkdir(parents=True, exist_ok=True)
        if os.path.isfile(out_file):
            out_dict = Props.read_from([out_file] + key_dict[key])
        else:
            out_dict = Props.read_from(key_dict[key])
        out_string = readable_json(out_dict)
        with open(out_file, "w") as w:
            w.write(out_string)

    os.system(f"rm {input_files}")


def build_file_dbs(input_files, output_dir):
    key_dict = get_keys(input_files)
    for key in list(key_dict):
        experiment, period, run, dtype = key.split("-")
        out_file = os.path.join(output_dir, f"{key}-filedb.h5")
        pathlib.Path(os.path.dirname(out_file)).mkdir(parents=True, exist_ok=True)
        file_path = f"{dtype}/{period}/{run}"
        cmd = f"{ut.runcmd(setup)} python3 -B {basedir}/scripts/build_fdb.py  --file_path {file_path} --output_file {out_file} --config {os.path.join(output_dir, 'file_db_config.json')}"
        os.system(cmd)


setup = snakemake.params.setup
basedir = snakemake.params.basedir

if os.getenv("PRODENV") in snakemake.params.filedb_path:
    file_db_config = {
        "data_dir": "$PRODENV",
        "tier_dirs": {
            "raw": ut.tier_raw_path(setup).replace(os.getenv("PRODENV"), ""),
            "dsp": ut.tier_dsp_path(setup).replace(os.getenv("PRODENV"), ""),
            "hit": ut.tier_hit_path(setup).replace(os.getenv("PRODENV"), ""),
            "pht": ut.tier_pht_path(setup).replace(os.getenv("PRODENV"), ""),
            "tcm": ut.tier_tcm_path(setup).replace(os.getenv("PRODENV"), ""),
            "evt": ut.tier_evt_path(setup).replace(os.getenv("PRODENV"), ""),
        },
        "file_format": {
            "raw": pat.get_pattern_tier(setup, "raw", check_in_cycle=False).replace(
                ut.tier_raw_path(setup), ""
            ),
            "dsp": pat.get_pattern_tier(setup, "dsp", check_in_cycle=False).replace(
                ut.tier_dsp_path(setup), ""
            ),
            "hit": pat.get_pattern_tier(setup, "hit", check_in_cycle=False).replace(
                ut.tier_hit_path(setup), ""
            ),
            "pht": pat.get_pattern_tier(setup, "pht", check_in_cycle=False).replace(
                ut.tier_hit_path(setup), ""
            ),
            "evt": pat.get_pattern_tier(setup, "evt", check_in_cycle=False).replace(
                ut.tier_evt_path(setup), ""
            ),
            "tcm": pat.get_pattern_tier(setup, "tcm", check_in_cycle=False).replace(
                ut.tier_tcm_path(setup), ""
            ),
        },
        "table_format": {
            "raw": "ch{ch:07d}/raw",
            "dsp": "ch{ch:07d}/dsp",
            "hit": "ch{ch:07d}/hit",
            "pht": "ch{ch:07d}/hit",
            "evt": "{grp}/evt",
            "tcm": "hardware_tcm_1",
        },
    }
else:
    file_db_config = {
        "data_dir": "/",
        "tier_dirs": {
            "raw": ut.tier_raw_path(setup),
            "dsp": ut.tier_dsp_path(setup),
            "hit": ut.tier_hit_path(setup),
            "tcm": ut.tier_tcm_path(setup),
            "pht": ut.tier_pht_path(setup),
            "evt": ut.tier_evt_path(setup),
        },
        "file_format": {
            "raw": pat.get_pattern_tier(setup, "raw", check_in_cycle=False).replace(
                ut.tier_raw_path(setup), ""
            ),
            "dsp": pat.get_pattern_tier(setup, "dsp", check_in_cycle=False).replace(
                ut.tier_dsp_path(setup), ""
            ),
            "hit": pat.get_pattern_tier(setup, "hit", check_in_cycle=False).replace(
                ut.tier_hit_path(setup), ""
            ),
            "pht": pat.get_pattern_tier(setup, "pht", check_in_cycle=False).replace(
                ut.tier_hit_path(setup), ""
            ),
            "evt": pat.get_pattern_tier(setup, "evt", check_in_cycle=False).replace(
                ut.tier_evt_path(setup), ""
            ),
            "tcm": pat.get_pattern_tier(setup, "tcm", check_in_cycle=False).replace(
                ut.tier_tcm_path(setup), ""
            ),
        },
        "table_format": {
            "raw": "ch{ch:07d}/raw",
            "dsp": "ch{ch:07d}/dsp",
            "hit": "ch{ch:07d}/hit",
            "pht": "ch{ch:07d}/hit",
            "evt": "{grp}/evt",
            "tcm": "hardware_tcm_1",
        },
    }

check_log_files(
    snakemake.params.log_path,
    snakemake.output.summary_log,
    snakemake.output.gen_output,
    warning_file=snakemake.output.warning_log,
)

if snakemake.wildcards.tier != "daq":
    os.makedirs(snakemake.params.filedb_path, exist_ok=True)
    with open(os.path.join(snakemake.params.filedb_path, "file_db_config.json"), "w") as w:
        json.dump(file_db_config, w, indent=2)

    build_file_dbs(snakemake.params.tmp_par_path, snakemake.params.filedb_path)
    os.system(f"rm {os.path.join(snakemake.params.filedb_path, 'file_db_config.json')}")

    build_valid_keys(snakemake.params.tmp_par_path, snakemake.params.valid_keys_path)

pathlib.Path(snakemake.output.gen_output).touch()
