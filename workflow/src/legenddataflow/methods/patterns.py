"""
This module contains all the patterns needed for the data production
"""

from __future__ import annotations

from pathlib import Path

from .paths import (
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
    return (
        "{experiment}-{period}-{run}-{datatype}-{timestamp}-{channel}-{processing_step}"
    )


def full_channel_pattern_with_extension():
    return "{experiment}-{period}-{run}-{datatype}-{timestamp}-{channel}-{processing_step}.{ext}"


def get_pattern_tier_daq_unsorted(setup, extension="orca"):
    if sandbox_path(setup) is not None:
        return Path(f"{sandbox_path(setup)}") / (
            "{experiment}-{period}-{run}-{datatype}-{timestamp}." + extension
        )
    return None


def get_pattern_tier_daq(setup, extension="orca", check_in_cycle=True):
    file_pattern = (
        Path(f"{tier_daq_path(setup)}")
        / "{datatype}"
        / "{period}"
        / "{run}"
        / ("{experiment}-{period}-{run}-{datatype}-{timestamp}." + extension)
    )
    if (
        tier_path(setup) not in str(file_pattern.resolve(strict=False))
        and check_in_cycle is True
    ):
        return "/tmp/" + file_pattern.name
    return file_pattern


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
            / (
                "{experiment}-{period}-{run}-{datatype}-{timestamp}-tier_"
                + f"{tier}.lh5"
            )
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
    if (
        tier_path(setup) not in str(file_pattern.resolve(strict=False))
        and check_in_cycle is True
    ):
        return "/tmp/" + file_pattern.name
    return file_pattern


def get_pattern_pars(
    setup, tier, name=None, datatype="cal", extension="yaml", check_in_cycle=True
):
    if datatype is None:
        datatype = "{datatype}"
    if tier in ["raw", "tcm", "dsp", "hit", "ann", "evt", "psp", "pht", "pan", "pet"]:
        if name is not None:
            file_pattern = (
                Path(get_pars_path(setup, tier))
                / datatype
                / "{period}"
                / "{run}"
                / (
                    "{experiment}-{period}-{run}-"
                    + datatype
                    + "-{timestamp}-par_"
                    + f"{tier}_{name}.{extension}"
                )
            )
        else:
            file_pattern = (
                Path(get_pars_path(setup, tier))
                / datatype
                / "{period}"
                / "{run}"
                / (
                    "{experiment}-{period}-{run}-"
                    + datatype
                    + "-{timestamp}-par_"
                    + f"{tier}.{extension}"
                )
            )
    else:
        msg = "invalid tier"
        raise Exception(msg)
    if (
        pars_path(setup) not in str(Path(file_pattern).resolve(strict=False))
        and check_in_cycle is True
    ):
        if name is None:
            return (
                "/tmp/{experiment}-{period}-{run}-"
                + datatype
                + "-{timestamp}-"
                + f"par_{tier}.{extension}"
            )
        return (
            "/tmp/{experiment}-{period}-{run}-" + datatype + "-{timestamp}-"
            f"par_{tier}_{name}.{extension}"
        )
    return file_pattern


def get_pattern_pars_inputs(setup, tier, name=None, ext="yaml"):
    if name is not None:
        return (
            Path(f"{par_overwrite_path(setup)}")
            / tier
            / "cal"
            / "{period}"
            / "{run}"
            / (
                "{experiment}-{period}-{run}-cal-{timestamp}-"
                + f"par_{tier}_{name}.{ext}"
            )
        )
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
            "{experiment}-{period}-{run}-"
            + datatype
            + "-{timestamp}-par_"
            + f"{tier}.{extension}"
        )
    return Path(f"{tmp_par_path(setup)}") / (
        "{experiment}-{period}-{run}-"
        + datatype
        + "-{timestamp}"
        + f"par_{tier}_{name}.{extension}"
    )


def get_pattern_pars_tmp_channel(
    setup, tier, name=None, datatype="cal", extension="yaml"
):
    if datatype is None:
        datatype = "{datatype}"
    if name is None:
        return Path(f"{tmp_par_path(setup)}") / (
            "{experiment}-{period}-{run}-"
            + datatype
            + "-{timestamp}-{channel}-par_"
            + f"{tier}.{extension}"
        )
    return Path(f"{tmp_par_path(setup)}") / (
        "{experiment}-{period}-{run}-"
        + datatype
        + "-{timestamp}-{channel}-par_"
        + f"{tier}_{name}.{extension}"
    )


def get_pattern_plts_tmp_channel(setup, tier, name=None, extension="pkl"):
    if name is None:
        return Path(f"{tmp_plts_path(setup)}") / (
            "{experiment}-{period}-{run}-cal-{timestamp}-{channel}-plt_"
            + f"{tier}.{extension}"
        )
    return Path(f"{tmp_plts_path(setup)}") / (
        "{experiment}-{period}-{run}-cal-{timestamp}-{channel}-plt_"
        + f"{tier}_{name}.{extension}"
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
    return (
        Path(f"{plts_path(setup)}")
        / tier
        / "cal"
        / "{period}"
        / "{run}"
        / (
            "{experiment}-{period}-{run}-cal-{timestamp}-plt_"
            + tier
            + "_"
            + name
            + ".dir"
        )
    )


def get_pattern_log(setup, processing_step, time):
    return (
        Path(f"{tmp_log_path(setup)}")
        / time
        / processing_step
        / (
            "{experiment}-{period}-{run}-{datatype}-{timestamp}-"
            + processing_step
            + ".log"
        )
    )


def get_pattern_log_channel(setup, processing_step, time, datatype="cal"):
    return (
        Path(f"{tmp_log_path(setup)}")
        / time
        / processing_step
        / (
            "{experiment}-{period}-{run}-"
            + datatype
            + "-{timestamp}-{channel}-"
            + processing_step
            + ".log"
        )
    )


def get_pattern_log_concat(setup, processing_step, time):
    return (
        Path(f"{tmp_log_path(setup)}")
        / time
        / processing_step
        / ("{experiment}-{period}-{run}-{datatype}-" + processing_step + ".log")
    )
