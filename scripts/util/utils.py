"""
This module contains all the utility needed for the data production.
These are mainly resolvers for the config.json dictionary,
and for substituting the pathvar within, also the conversion
from timestamp to unix time
"""

import copy
import os
import re
import shlex
import string
from datetime import datetime
from pathlib import Path

# from dateutil import parser

# For testing/debugging, use
# from scripts.utils import *
# import snakemake as smk
# setup = smk.load_configfile("config.json")["setups"]["l200"]


def sandbox_path(setup):
    if "sandbox_path" in setup["paths"]:
        return setup["paths"]["sandbox_path"]
    else:
        return None


def tier_daq_path(setup):
    return setup["paths"]["tier_daq"]


def tier_raw_blind_path(setup):
    return setup["paths"]["tier_raw_blind"]


def tier_path(setup):
    return setup["paths"]["tier"]


def tier_tcm_path(setup):
    return setup["paths"]["tier_tcm"]


def tier_raw_path(setup):
    return setup["paths"]["tier_raw"]


def tier_dsp_path(setup):
    return setup["paths"]["tier_dsp"]


def tier_hit_path(setup):
    return setup["paths"]["tier_hit"]


def tier_evt_path(setup):
    return setup["paths"]["tier_evt"]


def tier_psp_path(setup):
    return setup["paths"]["tier_psp"]


def tier_pht_path(setup):
    return setup["paths"]["tier_pht"]


def tier_pet_path(setup):
    return setup["paths"]["tier_pet"]


def tier_skm_path(setup):
    return setup["paths"]["tier_skm"]


def get_tier_path(setup, tier):
    if tier == "raw":
        return tier_raw_path(setup)
    elif tier == "tcm":
        return tier_tcm_path(setup)
    elif tier == "dsp":
        return tier_dsp_path(setup)
    elif tier == "hit":
        return tier_hit_path(setup)
    elif tier == "evt":
        return tier_evt_path(setup)
    elif tier == "psp":
        return tier_psp_path(setup)
    elif tier == "pht":
        return tier_pht_path(setup)
    elif tier == "pet":
        return tier_pet_path(setup)
    elif tier == "skm":
        return tier_skm_path(setup)
    else:
        msg = f"no tier matching:{tier}"
        raise ValueError(msg)


def config_path(setup):
    return setup["paths"]["config"]


def chan_map_path(setup):
    return setup["paths"]["chan_map"]


def metadata_path(setup):
    return setup["paths"]["metadata"]


def detector_db_path(setup):
    return setup["paths"]["detector_db"]


def par_raw_path(setup):
    return setup["paths"]["par_raw"]


def par_tcm_path(setup):
    return setup["paths"]["par_tcm"]


def par_dsp_path(setup):
    return setup["paths"]["par_dsp"]


def par_hit_path(setup):
    return setup["paths"]["par_hit"]


def par_evt_path(setup):
    return setup["paths"]["par_evt"]


def par_psp_path(setup):
    return setup["paths"]["par_psp"]


def par_pht_path(setup):
    return setup["paths"]["par_pht"]


def par_pet_path(setup):
    return setup["paths"]["par_pet"]


def pars_path(setup):
    return setup["paths"]["par"]


def get_pars_path(setup, tier):
    if tier == "raw":
        return par_raw_path(setup)
    elif tier == "tcm":
        return par_tcm_path(setup)
    elif tier == "dsp":
        return par_dsp_path(setup)
    elif tier == "hit":
        return par_hit_path(setup)
    elif tier == "evt":
        return par_evt_path(setup)
    elif tier == "psp":
        return par_psp_path(setup)
    elif tier == "pht":
        return par_pht_path(setup)
    elif tier == "pet":
        return par_pet_path(setup)
    else:
        msg = f"no tier matching:{tier}"
        raise ValueError(msg)


def tmp_par_path(setup):
    return setup["paths"]["tmp_par"]


def tmp_plts_path(setup):
    return setup["paths"]["tmp_plt"]


def plts_path(setup):
    return setup["paths"]["plt"]


def par_overwrite_path(setup):
    return setup["paths"]["par_overwrite"]


def log_path(setup):
    return setup["paths"]["log"]


def tmp_log_path(setup):
    return setup["paths"]["tmp_log"]


def filelist_path(setup):
    return setup["paths"]["tmp_filelists"]


def runcmd(setup, aslist=False):
    cmdline = shlex.split(setup["execenv"]["cmd"])
    cmdline += ["--env=" + "'PYTHONUSERBASE=" + f"{setup['paths']['install']}" + "'"]
    if "env" in setup["execenv"]:
        cmdline += [f'--env="{var}={val}"' for var, val in setup["execenv"]["env"].items()]

    cmdline += shlex.split(setup["execenv"]["arg"])

    if aslist:
        return cmdline

    return " ".join(cmdline)


def subst_vars_impl(x, var_values, ignore_missing=False):
    if isinstance(x, str):
        if "$" in x:
            if ignore_missing:
                return string.Template(x).safe_substitute(var_values)
            else:
                return string.Template(x).substitute(var_values)
        else:
            return x
    if isinstance(x, dict):
        for key in x:
            value = x[key]
            new_value = subst_vars_impl(value, var_values, ignore_missing)
            if new_value is not value:
                x[key] = new_value
        return x
    if isinstance(x, list):
        for i in range(len(x)):
            value = x[i]
            new_value = subst_vars_impl(value, var_values, ignore_missing)
            if new_value is not value:
                x[i] = new_value
        return x
    else:
        return x


def subst_vars(props, var_values=None, use_env=False, ignore_missing=False):
    if var_values is None:
        var_values = {}
    combined_var_values = var_values
    if use_env:
        combined_var_values = dict(iter(os.environ.items()))
        combined_var_values.update(copy.copy(var_values))
    subst_vars_impl(props, combined_var_values, ignore_missing)


def subst_vars_in_snakemake_config(workflow, config):
    config_filename = workflow.overwrite_configfiles[0]  # ToDo: Better way of handling this?
    subst_vars(
        config,
        var_values={"_": os.path.dirname(config_filename)},
        use_env=True,
        ignore_missing=False,
    )


def run_splitter(files):
    """
    Returns list containing lists of each run
    """

    runs = []
    run_files = []
    for file in files:
        base = os.path.basename(file)
        file_name = os.path.splitext(base)[0]
        parts = file_name.split("-")
        run_no = parts[3]
        if run_no not in runs:
            runs.append(run_no)
            run_files.append([])
        for i, run in enumerate(runs):
            if run == run_no:
                run_files[i].append(file)
    return run_files


def unix_time(value):
    if isinstance(value, str):
        return datetime.timestamp(datetime.strptime(value, "%Y%m%dT%H%M%SZ"))
    else:
        msg = f"Can't convert type {type(value)} to unix time"
        raise ValueError(msg)


def set_last_rule_name(workflow, new_name):
    """Sets the name of the most recently created rule to be `new_name`.
    Useful when creating rules dynamically (i.e. unnamed).

    Warning
    -------
    This could mess up the workflow. Use at your own risk.
    """
    rules = workflow._rules
    last_key = next(reversed(rules))
    assert last_key == rules[last_key].name

    rules[new_name] = rules.pop(last_key)
    rules[new_name].name = new_name

    if workflow.default_target == last_key:
        workflow.default_target = new_name

    if last_key in workflow._localrules:
        workflow._localrules.remove(last_key)
        workflow._localrules.add(new_name)

    workflow.check_localrules()


def as_ro(config, path):
    if "read_only_fs_sub_pattern" not in config or config["read_only_fs_sub_pattern"] is None:
        return path

    sub_pattern = config["read_only_fs_sub_pattern"]

    if isinstance(path, str):
        return re.sub(*sub_pattern, path)
    if isinstance(path, Path):
        return Path(re.sub(*sub_pattern, path.name))

    return [as_ro(config, p) for p in path]
