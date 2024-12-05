"""
This module contains all the patterns needed for the data production
"""

from pathlib import Path

from .utils import (
    get_pars_path,
    get_tier_path,
    par_overwrite_path,
    pars_path,
    plts_path,
    sandbox_path,
    tier_daq_path,
    tier_path,
    tier_raw_blind_path,
    tmp_log_path,
    tmp_par_path,
    tmp_plts_path,
)


# key_mask
def key_pattern():
    return "{experiment}-{period}-{run}-{datatype}-{timestamp}"


def processing_pattern():
    return key_pattern() + "-{processing_step}.{ext}"


def par_validity_pattern():
    return "{datatype}/{period}/{run}/" + processing_pattern()


def par_pattern():
    return "{experiment}-{period}-{run}-{datatype}-{timestamp}-par_{name}"


def par_overwrite_pattern():
    return key_pattern() + "-{processing_step}-overwrite"


def processing_overwrite_pattern():
    return par_overwrite_pattern() + ".{ext}"


def full_channel_pattern():
    return "{experiment}-{period}-{run}-{datatype}-{timestamp}-{channel}-{processing_step}"


def full_channel_pattern_with_extension():
    return "{experiment}-{period}-{run}-{datatype}-{timestamp}-{channel}-{processing_step}.{ext}"


def get_pattern_unsorted_data(setup):
    if sandbox_path(setup) is not None:
        return (
            Path(f"{sandbox_path(setup)}")
            / "{experiment}-{period}-{run}-{datatype}-{timestamp}.orca"
        )
    else:
        return None


def get_pattern_tier_daq(setup):
    return (
        Path(f"{tier_daq_path(setup)}")
        / "{datatype}"
        / "{period}"
        / "{run}"
        / "{experiment}-{period}-{run}-{datatype}-{timestamp}.orca"
    )


def get_pattern_tier_raw_blind(setup):
    return (
        Path(f"{tier_raw_blind_path(setup)}")
        / "phy"
        / "{period}"
        / "{run}"
        / "{experiment}-{period}-{run}-phy-{timestamp}-tier_raw.lh5"
    )


def get_pattern_tier(setup, tier, check_in_cycle=True):
    if tier in ["raw", "tcm", "dsp", "hit", "ann", "evt", "psp", "pht", "pan", "pet"]:
        file_pattern = (
            Path(get_tier_path(setup, tier))
            / "{datatype}"
            / "{period}"
            / "{run}"
            / ("{experiment}-{period}-{run}-{datatype}-{timestamp}-tier_" + f"{tier}.lh5")
        )
    elif tier in ["evt_concat", "pet_concat"]:
        file_pattern = (
            Path(get_tier_path(setup, tier[:3]))
            / "{datatype}"
            / ("{experiment}-{period}-{run}-{datatype}-tier_" + f"{tier[:3]}.lh5")
        )

    elif tier == "skm":
        file_pattern = (
            Path(f"{get_tier_path(setup, tier)}")
            / "phy"
            / "{experiment}-{period}-{run}-{datatype}-tier_skm.lh5"
        )
    else:
        msg = "invalid tier"
        raise Exception(msg)
    if tier_path(setup) not in str(file_pattern.resolve(strict=False)) and check_in_cycle is True:
        return "/tmp/" + file_pattern.name
    else:
        return file_pattern


def get_pattern_pars(setup, tier, name=None, extension="yaml", check_in_cycle=True):
    if tier in ["raw", "tcm", "dsp", "hit", "ann", "evt", "psp", "pht", "pan", "pet"]:
        if name is not None:
            return (
                Path(get_pars_path(setup, tier))
                / "cal"
                / "{period}"
                / "{run}"
                / (
                    "{experiment}-{period}-{run}-cal-{timestamp}-par_"
                    + f"{tier}_{name}.{extension}"
                )
            )
        else:
            file_pattern = (
                Path(get_pars_path(setup, tier))
                / "cal"
                / "{period}"
                / "{run}"
                / ("{experiment}-{period}-{run}-cal-{timestamp}-par_" + f"{tier}.{extension}")
            )
    else:
        msg = "invalid tier"
        raise Exception(msg)
    if (
        pars_path(setup) not in str(Path(file_pattern).resolve(strict=False))
        and check_in_cycle is True
    ):
        if name is None:
            return "/tmp/{experiment}-{period}-{run}-cal-{timestamp}-" + f"par_{tier}.{extension}"
        else:
            return (
                "/tmp/{experiment}-{period}-{run}-cal-{timestamp}-"
                f"par_{tier}_{name}.{extension}"
            )
    else:
        return file_pattern


def get_pattern_pars_inputs(setup, tier, name=None, ext="yaml"):
    if name is not None:
        return (
            Path(f"{par_overwrite_path(setup)}")
            / tier
            / "cal"
            / "{period}"
            / "{run}"
            / ("{experiment}-{period}-{run}-cal-{timestamp}-" + f"par_{tier}_{name}.{ext}")
        )
    else:
        return (
            Path(f"{par_overwrite_path(setup)}")
            / tier
            / "cal"
            / "{period}"
            / "{run}"
            / ("{experiment}-{period}-{run}-cal-{timestamp}-" + f"par_{tier}.{ext}")
        )


def get_pattern_pars_overwrite(setup, tier, name=None, extension="yaml"):
    if name is not None:
        return (
            Path(f"{par_overwrite_path(setup)}")
            / tier
            / "{datatype}"
            / "{period}"
            / "{run}"
            / (
                "{experiment}-{period}-{run}-{datatype}-{timestamp}-par_"
                f"{tier}_{name}-overwrite.{extension}"
            )
        )
    else:
        return (
            Path(f"{par_overwrite_path(setup)}")
            / tier
            / "{datatype}"
            / "{period}"
            / "{run}"
            / (
                "{experiment}-{period}-{run}-{datatype}-{timestamp}-par_"
                + tier
                + f"-overwrite.{extension}"
            )
        )


def get_pattern_pars_tmp(setup, tier, name=None, datatype=None, extension="yaml"):
    if datatype is None:
        datatype = "{datatype}"
    if name is None:
        return Path(f"{tmp_par_path(setup)}") / (
            "{experiment}-{period}-{run}-" + datatype + "-{timestamp}-par_" + f"{tier}.{extension}"
        )
    else:
        return Path(f"{tmp_par_path(setup)}") / (
            "{experiment}-{period}-{run}-"
            + datatype
            + "-{timestamp}"
            + f"par_{tier}_{name}.{extension}"
        )


def get_pattern_pars_tmp_channel(setup, tier, name=None, extension="yaml"):
    if name is None:
        return Path(f"{tmp_par_path(setup)}") / (
            "{experiment}-{period}-{run}-cal-{timestamp}-{channel}-par_" + f"{tier}.{extension}"
        )
    else:
        return Path(f"{tmp_par_path(setup)}") / (
            "{experiment}-{period}-{run}-cal-{timestamp}-{channel}-par_"
            + f"{tier}_{name}.{extension}"
        )


def get_pattern_plts_tmp_channel(setup, tier, name=None):
    if name is None:
        return Path(f"{tmp_plts_path(setup)}") / (
            "{experiment}-{period}-{run}-cal-{timestamp}-{channel}-plt_" + tier + ".pkl"
        )
    else:
        return Path(f"{tmp_plts_path(setup)}") / (
            "{experiment}-{period}-{run}-cal-{timestamp}-{channel}-plt_" + f"{tier}_{name}.pkl"
        )


def get_pattern_plts(setup, tier, name=None):
    if name is None:
        return (
            Path(f"{plts_path(setup)}")
            / tier
            / "cal"
            / "{period}"
            / "{run}"
            / ("{experiment}-{period}-{run}-cal-{timestamp}-plt_" + tier + ".dir")
        )
    else:
        return (
            Path(f"{plts_path(setup)}")
            / tier
            / "cal"
            / "{period}"
            / "{run}"
            / ("{experiment}-{period}-{run}-cal-{timestamp}-plt_" + tier + "_" + name + ".dir")
        )


def get_pattern_log(setup, processing_step):
    return (
        Path(f"{tmp_log_path(setup)}")
        / processing_step
        / ("{experiment}-{period}-{run}-{datatype}-{timestamp}-" + processing_step + ".log")
    )


def get_pattern_log_channel(setup, processing_step):
    return (
        Path(f"{tmp_log_path(setup)}")
        / processing_step
        / ("{experiment}-{period}-{run}-cal-{timestamp}-{channel}-" + processing_step + ".log")
    )


def get_pattern_log_concat(setup, processing_step):
    return (
        Path(f"{tmp_log_path(setup)}")
        / processing_step
        / ("{experiment}-{period}-{run}-{datatype}-" + processing_step + ".log")
    )
