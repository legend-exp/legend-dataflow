"""
This module contains resolvers for the config.json dictionary
"""

from __future__ import annotations


def sandbox_path(setup):
    if "sandbox_path" in setup["paths"]:
        return setup["paths"]["sandbox_path"]
    return None


def tier_daq_path(setup):
    return setup["paths"]["tier_daq"]


def tier_raw_blind_path(setup):
    return setup["paths"]["tier_raw_blind"]


def tier_path(setup):
    return setup["paths"]["tier"]


def get_tier_path(setup, tier):
    if tier in [
        "raw",
        "tcm",
        "dsp",
        "hit",
        "ann",
        "evt",
        "psp",
        "pht",
        "pan",
        "pet",
        "skm",
    ]:
        return setup["paths"][f"tier_{tier}"]
    msg = f"no tier matching:{tier}"
    raise ValueError(msg)


def pars_path(setup):
    return setup["paths"]["par"]


def get_pars_path(setup, tier):
    if tier in ["raw", "tcm", "dsp", "hit", "evt", "psp", "pht", "pet"]:
        return setup["paths"][f"par_{tier}"]
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


def config_path(setup):
    return setup["paths"]["config"]


def chan_map_path(setup):
    return setup["paths"]["chan_map"]


def det_status_path(setup):
    return setup["paths"]["detector_status"]


def metadata_path(setup):
    return setup["paths"]["metadata"]


def detector_db_path(setup):
    return setup["paths"]["detector_db"]


def log_path(setup):
    return setup["paths"]["log"]


def tmp_log_path(setup):
    return setup["paths"]["tmp_log"]


def filelist_path(setup):
    return setup["paths"]["tmp_filelists"]
