"""
Helper functions for running data production
"""

import pathlib, os
from scripts.util.patterns import (
    par_overwrite_path,
    par_raw_path,
    get_pattern_unsorted_data,
    get_pattern_tier_daq,
    get_pattern_tier_raw,
    get_pattern_plts_tmp_channel,
)
from scripts.util import ProcessingFileKey
from scripts.util import utils


def ro(path):
    return utils.as_ro(setup, path)


def get_blinding_curve_file(wildcards):
    """func to get the blinding calibration curves from the overrides"""
    par_files = pars_catalog.get_calib_files(
        Path(par_overwrite_path(setup)) / "raw" / "validity.jsonl",
        wildcards.timestamp,
    )
    if isinstance(par_files, str):
        return str(Path(par_overwrite_path(setup)) / "raw" / par_files)
    else:
        return [
            str(Path(par_overwrite_path(setup)) / "raw" / par_file)
            for par_file in par_files
        ]


def get_blinding_check_file(wildcards):
    """func to get the right blinding check file"""
    par_files = pars_catalog.get_calib_files(
        Path(par_raw_path(setup)) / "validity.jsonl", wildcards.timestamp
    )
    if isinstance(par_files, str):
        return str(Path(par_raw_path(setup)) / par_files)
    else:
        return [str(Path(par_raw_path(setup)) / par_file) for par_file in par_files]


# def get_pattern(tier):
#     """
#     This func gets the search pattern for the relevant tier passed.
#     """
#     if tier == "daq":
#         return get_pattern_unsorted_data(setup)
#     elif tier == "raw":
#         return get_pattern_tier_daq(setup)
#     else:
#         return get_pattern_tier_raw(setup)


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


def get_svm_file(wildcards, tier, name):
    par_overwrite_file = os.path.join(par_overwrite_path(setup), tier, "validity.jsonl")
    pars_files_overwrite = pars_catalog.get_calib_files(
        par_overwrite_file, wildcards.timestamp
    )
    for pars_file in pars_files_overwrite:
        if name in pars_file:
            return os.path.join(par_overwrite_path(setup), tier, pars_file)
    raise ValueError(f"Could not find model in {pars_files_overwrite}")


def get_overwrite_file(tier, wildcards=None, timestamp=None, name=None):
    par_overwrite_file = os.path.join(par_overwrite_path(setup), tier, "validity.jsonl")
    if timestamp is not None:
        pars_files_overwrite = pars_catalog.get_calib_files(
            par_overwrite_file, timestamp
        )
    else:
        pars_files_overwrite = pars_catalog.get_calib_files(
            par_overwrite_file, wildcards.timestamp
        )
    if name is None:
        fullname = f"{tier}-overwrite.json"
    else:
        fullname = f"{tier}_{name}-overwrite.json"
    out_files = []
    for pars_file in pars_files_overwrite:
        if fullname in pars_file:
            out_files.append(os.path.join(par_overwrite_path(setup), tier, pars_file))
    if len(out_files) == 0:
        raise ValueError(f"Could not find name in {pars_files_overwrite}")
    else:
        return out_files
