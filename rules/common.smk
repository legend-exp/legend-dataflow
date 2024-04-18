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


def read_filelist(wildcards):
    with checkpoints.gen_filelist.get(
        label=wildcards.label, tier=wildcards.tier, extension="file"
    ).output[0].open() as f:
        files = f.read().splitlines()
        return files


def read_filelist_phy(wildcards, tier):
    label = f"all-{wildcards.experiment}-{wildcards.period}-{wildcards.run}-phy"
    with checkpoints.gen_filelist.get(label=label, tier=tier, extension="file").output[
        0
    ].open() as f:
        files = f.read().splitlines()
        return files


def read_filelist_cal(wildcards, tier):
    label = f"all-{wildcards.experiment}-{wildcards.period}-{wildcards.run}-cal"
    with checkpoints.gen_filelist.get(label=label, tier=tier, extension="file").output[
        0
    ].open() as f:
        files = f.read().splitlines()
        return files


def read_filelist_fft(wildcards, tier):
    label = f"all-{wildcards.experiment}-{wildcards.period}-{wildcards.run}-fft"
    with checkpoints.gen_filelist.get(label=label, tier=tier, extension="file").output[
        0
    ].open() as f:
        files = f.read().splitlines()
        return files


def read_filelist_pars_cal_channel(wildcards, tier):
    """
    This function will read the filelist of the channels and return a list of dsp files one for each channel
    """
    label = f"all-{wildcards.experiment}-{wildcards.period}-{wildcards.run}-cal-{wildcards.timestamp}-channels"
    with checkpoints.gen_filelist.get(
        label=label, tier=f"par_{tier}", extension="chan"
    ).output[0].open() as f:
        files = f.read().splitlines()
        return files


def read_filelist_plts_cal_channel(wildcards, tier):
    """
    This function will read the filelist of the channels and return a list of dsp files one for each channel
    """
    label = f"all-{wildcards.experiment}-{wildcards.period}-{wildcards.run}-cal-{wildcards.timestamp}-channels"
    with checkpoints.gen_filelist.get(
        label=label, tier=f"plt_{tier}", extension="chan"
    ).output[0].open() as f:
        files = f.read().splitlines()
        return files


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


def get_pattern(tier):
    """
    This func gets the search pattern for the relevant tier passed.
    """
    if tier == "daq":
        return get_pattern_unsorted_data(setup)
    elif tier == "raw":
        return get_pattern_tier_daq(setup)
    else:
        return get_pattern_tier_raw(setup)


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
