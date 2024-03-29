"""
This module contains all the patterns needed for the data production
"""

import os

from .utils import (
    par_dsp_path,
    par_evt_path,
    par_hit_path,
    par_overwrite_path,
    par_pht_path,
    par_psp_path,
    par_raw_path,
    par_tcm_path,
    pars_path,
    plts_path,
    sandbox_path,
    tier_daq_path,
    tier_dsp_path,
    tier_evt_path,
    tier_hit_path,
    tier_path,
    tier_pet_path,
    tier_pht_path,
    tier_psp_path,
    tier_raw_blind_path,
    tier_raw_path,
    tier_skm_path,
    tier_tcm_path,
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
        return os.path.join(
            f"{sandbox_path(setup)}",
            "{experiment}-{period}-{run}-{datatype}-{timestamp}.orca",
        )
    else:
        return None


def get_pattern_tier_daq(setup):
    return os.path.join(
        f"{tier_daq_path(setup)}",
        "{datatype}",
        "{period}",
        "{run}",
        "{experiment}-{period}-{run}-{datatype}-{timestamp}.orca",
    )


def get_pattern_tier_raw(setup):
    return os.path.join(
        f"{tier_raw_path(setup)}",
        "{datatype}",
        "{period}",
        "{run}",
        "{experiment}-{period}-{run}-{datatype}-{timestamp}-tier_raw.lh5",
    )


def get_pattern_tier_raw_blind(setup):
    return os.path.join(
        f"{tier_raw_blind_path(setup)}",
        "phy",
        "{period}",
        "{run}",
        "{experiment}-{period}-{run}-phy-{timestamp}-tier_raw.lh5",
    )


def get_pattern_tier_tcm(setup):
    return os.path.join(
        f"{tier_tcm_path(setup)}",
        "{datatype}",
        "{period}",
        "{run}",
        "{experiment}-{period}-{run}-{datatype}-{timestamp}-tier_tcm.lh5",
    )


def get_pattern_tier_dsp(setup):
    return os.path.join(
        f"{tier_dsp_path(setup)}",
        "{datatype}",
        "{period}",
        "{run}",
        "{experiment}-{period}-{run}-{datatype}-{timestamp}-tier_dsp.lh5",
    )


def get_pattern_tier_hit(setup):
    return os.path.join(
        f"{tier_hit_path(setup)}",
        "{datatype}",
        "{period}",
        "{run}",
        "{experiment}-{period}-{run}-{datatype}-{timestamp}-tier_hit.lh5",
    )


def get_pattern_tier_evt(setup):
    return os.path.join(
        f"{tier_evt_path(setup)}",
        "{datatype}",
        "{period}",
        "{run}",
        "{experiment}-{period}-{run}-{datatype}-{timestamp}-tier_evt.lh5",
    )


def get_pattern_tier_psp(setup):
    return os.path.join(
        f"{tier_psp_path(setup)}",
        "{datatype}",
        "{period}",
        "{run}",
        "{experiment}-{period}-{run}-{datatype}-{timestamp}-tier_psp.lh5",
    )


def get_pattern_tier_pht(setup):
    return os.path.join(
        f"{tier_pht_path(setup)}",
        "{datatype}",
        "{period}",
        "{run}",
        "{experiment}-{period}-{run}-{datatype}-{timestamp}-tier_pht.lh5",
    )


def get_pattern_tier_pet(setup):
    return os.path.join(
        f"{tier_pet_path(setup)}",
        "{datatype}",
        "{period}",
        "{run}",
        "{experiment}-{period}-{run}-{datatype}-{timestamp}-tier_pet.lh5",
    )


def get_pattern_tier_skm(setup):
    return os.path.join(
        f"{tier_skm_path(setup)}",
        "{datatype}",
        "{period}",
        "{run}",
        "{experiment}-{period}-{run}-{datatype}-{timestamp}-tier_skm.lh5",
    )


def get_pattern_tier(setup, tier, check_in_cycle=True):
    if tier == "daq":
        file_pattern = get_pattern_tier_daq(setup)
    elif tier == "raw":
        file_pattern = get_pattern_tier_raw(setup)
    elif tier == "tcm":
        file_pattern = get_pattern_tier_tcm(setup)
    elif tier == "dsp":
        file_pattern = get_pattern_tier_dsp(setup)
    elif tier == "hit":
        file_pattern = get_pattern_tier_hit(setup)
    elif tier == "evt":
        file_pattern = get_pattern_tier_evt(setup)
    elif tier == "psp":
        file_pattern = get_pattern_tier_psp(setup)
    elif tier == "pht":
        file_pattern = get_pattern_tier_pht(setup)
    elif tier == "pet":
        file_pattern = get_pattern_tier_pet(setup)
    elif tier == "skm":
        file_pattern = get_pattern_tier_skm(setup)
    else:
        msg = "invalid tier"
        raise Exception(msg)
    if (
        tier_path(setup) not in file_pattern
        and check_in_cycle is True
        and ".." not in file_pattern
    ):
        return "/tmp/{experiment}-{period}-{run}-{datatype}-{timestamp}" + f"tier_{tier}.lh5"
    else:
        return file_pattern


def get_pattern_par_raw(setup, name=None, extension="json"):
    if name is not None:
        return os.path.join(
            f"{par_raw_path(setup)}",
            "cal",
            "{period}",
            "{run}",
            "{experiment}-{period}-{run}-cal-{timestamp}-par_raw_" + f"{name}.{extension}",
        )
    else:
        return os.path.join(
            f"{par_raw_path(setup)}",
            "cal",
            "{period}",
            "{run}",
            "{experiment}-{period}-{run}-cal-{timestamp}-par_raw" + f".{extension}",
        )


def get_pattern_par_tcm(setup, name=None, extension="json"):
    if name is not None:
        return os.path.join(
            f"{par_tcm_path(setup)}",
            "cal",
            "{period}",
            "{run}",
            "{experiment}-{period}-{run}-cal-{timestamp}-par_tcm_" + f"{name}.{extension}",
        )
    else:
        return os.path.join(
            f"{par_tcm_path(setup)}",
            "cal",
            "{period}",
            "{run}",
            "{experiment}-{period}-{run}-cal-{timestamp}-par_tcm" + f".{extension}",
        )


def get_pattern_par_dsp(setup, name=None, extension="json"):
    if name is not None:
        return os.path.join(
            f"{par_dsp_path(setup)}",
            "cal",
            "{period}",
            "{run}",
            "{experiment}-{period}-{run}-cal-{timestamp}-par_dsp_" + f"{name}.{extension}",
        )
    else:
        return os.path.join(
            f"{par_dsp_path(setup)}",
            "cal",
            "{period}",
            "{run}",
            "{experiment}-{period}-{run}-cal-{timestamp}-par_dsp" + f".{extension}",
        )


def get_pattern_par_hit(setup, name=None, extension="json"):
    if name is not None:
        return os.path.join(
            f"{par_hit_path(setup)}",
            "cal",
            "{period}",
            "{run}",
            "{experiment}-{period}-{run}-cal-{timestamp}-par_hit_" + f"{name}.{extension}",
        )
    else:
        return os.path.join(
            f"{par_hit_path(setup)}",
            "cal",
            "{period}",
            "{run}",
            "{experiment}-{period}-{run}-cal-{timestamp}-par_hit" + f".{extension}",
        )


def get_pattern_par_evt(setup, name=None, extension="json"):
    if name is not None:
        return os.path.join(
            f"{par_evt_path(setup)}",
            "cal",
            "{period}",
            "{run}",
            "{experiment}-{period}-{run}-cal-{timestamp}-par_evt_" + f"{name}.{extension}",
        )
    else:
        return os.path.join(
            f"{par_evt_path(setup)}",
            "cal",
            "{period}",
            "{run}",
            "{experiment}-{period}-{run}-cal-{timestamp}-par_evt" + f".{extension}",
        )


def get_pattern_par_psp(setup, name=None, extension="json"):
    if name is not None:
        return os.path.join(
            f"{par_psp_path(setup)}",
            "cal",
            "{period}",
            "{run}",
            "{experiment}-{period}-{run}-cal-{timestamp}-par_psp_" + f"{name}.{extension}",
        )
    else:
        return os.path.join(
            f"{par_psp_path(setup)}",
            "cal",
            "{period}",
            "{run}",
            "{experiment}-{period}-{run}-cal-{timestamp}-par_psp" + f".{extension}",
        )


def get_pattern_par_pht(setup, name=None, extension="json"):
    if name is not None:
        return os.path.join(
            f"{par_pht_path(setup)}",
            "cal",
            "{period}",
            "{run}",
            "{experiment}-{period}-{run}-cal-{timestamp}-par_pht_" + f"{name}.{extension}",
        )
    else:
        return os.path.join(
            f"{par_pht_path(setup)}",
            "cal",
            "{period}",
            "{run}",
            "{experiment}-{period}-{run}-cal-{timestamp}-par_pht" + f".{extension}",
        )


def get_pattern_par_pet(setup, name=None, extension="json"):
    if name is not None:
        return os.path.join(
            f"{par_evt_path(setup)}",
            "cal",
            "{period}",
            "{run}",
            "{experiment}-{period}-{run}-cal-{timestamp}-par_pet_" + f"{name}.{extension}",
        )
    else:
        return os.path.join(
            f"{par_evt_path(setup)}",
            "cal",
            "{period}",
            "{run}",
            "{experiment}-{period}-{run}-cal-{timestamp}-par_pet" + f".{extension}",
        )


def get_pattern_pars(setup, tier, name=None, extension="json", check_in_cycle=True):
    if tier == "raw":
        file_pattern = get_pattern_par_raw(setup, name, extension)
    elif tier == "tcm":
        file_pattern = get_pattern_par_tcm(setup, name, extension)
    elif tier == "dsp":
        file_pattern = get_pattern_par_dsp(setup, name, extension)
    elif tier == "hit":
        file_pattern = get_pattern_par_hit(setup, name, extension)
    elif tier == "evt":
        file_pattern = get_pattern_par_evt(setup, name, extension)
    elif tier == "psp":
        file_pattern = get_pattern_par_psp(setup, name, extension)
    elif tier == "pht":
        file_pattern = get_pattern_par_pht(setup, name, extension)
    elif tier == "pet":
        file_pattern = get_pattern_par_pet(setup, name, extension)
    else:
        msg = "invalid tier"
        raise Exception(msg)
    if (
        pars_path(setup) not in file_pattern
        and check_in_cycle is True
        and ".." not in file_pattern
    ):
        if name is None:
            return "/tmp/{experiment}-{period}-{run}-cal-{timestamp}" + f"par_{tier}.{extension}"
        else:
            return (
                "/tmp/{experiment}-{period}-{run}-cal-{timestamp}"
                + f"par_{tier}_{name}.{extension}"
            )
    else:
        return file_pattern


def get_pattern_pars_overwrite(setup, tier, name=None):
    if name is not None:
        return os.path.join(
            f"{par_overwrite_path(setup)}",
            tier,
            "{datatype}",
            "{period}",
            "{run}",
            "{experiment}-{period}-{run}-{datatype}-{timestamp}-par_"
            + tier
            + "_"
            + name
            + "-overwrite.json",
        )
    else:
        return os.path.join(
            f"{par_overwrite_path(setup)}",
            tier,
            "{datatype}",
            "{period}",
            "{run}",
            "{experiment}-{period}-{run}-{datatype}-{timestamp}-par_" + tier + "-overwrite.json",
        )


def get_pattern_pars_tmp(setup, tier, name=None, datatype=None):
    if datatype is None:
        datatype = "{datatype}"
    if name is None:
        return os.path.join(
            f"{tmp_par_path(setup)}",
            "{experiment}-{period}-{run}-" + datatype + "-{timestamp}-par_" + tier + ".json",
        )
    else:
        return os.path.join(
            f"{tmp_par_path(setup)}",
            "{experiment}-{period}-{run}-"
            + datatype
            + "-{timestamp}-par_"
            + tier
            + "_"
            + name
            + ".json",
        )


def get_pattern_pars_tmp_channel(setup, tier, name=None, extension="json"):
    if name is None:
        return os.path.join(
            f"{tmp_par_path(setup)}",
            "{experiment}-{period}-{run}-cal-{timestamp}-{channel}-par_" + f"{tier}.{extension}",
        )
    else:
        return os.path.join(
            f"{tmp_par_path(setup)}",
            "{experiment}-{period}-{run}-cal-{timestamp}-{channel}-par_"
            + f"{tier}_{name}.{extension}",
        )


def get_pattern_plts_tmp_channel(setup, tier, name=None):
    if name is None:
        return os.path.join(
            f"{tmp_plts_path(setup)}",
            "{experiment}-{period}-{run}-cal-{timestamp}-{channel}-plt_" + tier + ".pkl",
        )
    else:
        return os.path.join(
            f"{tmp_plts_path(setup)}",
            "{experiment}-{period}-{run}-cal-{timestamp}-{channel}-plt_"
            + tier
            + "_"
            + name
            + ".pkl",
        )


def get_pattern_plts(setup, tier, name=None):
    if name is None:
        return os.path.join(
            f"{plts_path(setup)}",
            tier,
            "cal",
            "{period}",
            "{run}",
            "{experiment}-{period}-{run}-cal-{timestamp}-plt_" + tier + ".dir",
        )
    else:
        return os.path.join(
            f"{plts_path(setup)}",
            tier,
            "cal",
            "{period}",
            "{run}",
            "{experiment}-{period}-{run}-cal-{timestamp}-plt_" + tier + "_" + name + ".dir",
        )


def get_energy_grids_pattern_combine(setup):
    return os.path.join(
        f"{tmp_par_path(setup)}",
        "dsp",
        "cal",
        "{{period}}",
        "{{run}}",
        "par_dsp_energy_grid",
        "{{channel}}",
        "{{experiment}}-{{period}}-{{run}}-cal-{{timestamp}}-{{channel}}-{peak}-par_dsp_energy_grid.pkl",
    )


def get_pattern_log(setup, processing_step):
    return os.path.join(
        f"{tmp_log_path(setup)}",
        processing_step,
        "{experiment}-{period}-{run}-{datatype}-{timestamp}-" + processing_step + ".log",
    )


def get_pattern_log_channel(setup, processing_step):
    return os.path.join(
        f"{tmp_log_path(setup)}",
        processing_step,
        "{experiment}-{period}-{run}-cal-{timestamp}-{channel}-" + processing_step + ".log",
    )
