import os
from .utils import *

def key_pattern():
    return "{experiment}-{period}-{run}-{datatype}-{timestamp}"

def processing_pattern():
    return key_pattern()+'-{processing_step}'

def par_pattern():
    return "{experiment}-{period}-{run}-{datatype}-{timestamp}-par_{name}"

def par_overwrite_pattern():
    return par_pattern()+"-overwrite"

def processing_overwrite_pattern():
    return processing_pattern()+"-overwrite"

def full_channel_pattern():
    return "{experiment}-{period}-{run}-{datatype}-{timestamp}-{channel}-{processing_step}"

def get_pattern_tier_daq(setup):
    return os.path.join(f"{tier_daq_path(setup)}", "{datatype}", "{period}", "{run}", "{experiment}-{period}-{run}-{datatype}-{timestamp}.orca")

def get_pattern_tier_raw(setup):
    return os.path.join(f"{tier_raw_path(setup)}", "{datatype}","{period}", "{run}", "{experiment}-{period}-{run}-{datatype}-{timestamp}-tier_raw.lh5")

def get_pattern_tier_tcm(setup):
    return os.path.join(f"{tier_tcm_path(setup)}", "{datatype}","{period}", "{run}", "{experiment}-{period}-{run}-{datatype}-{timestamp}-tier_tcm.lh5")

def get_pattern_tier_dsp(setup):
    return os.path.join(f"{tier_dsp_path(setup)}", "{datatype}","{period}", "{run}", "{experiment}-{period}-{run}-{datatype}-{timestamp}-tier_dsp.lh5")

def get_pattern_tier_hit(setup):
    return os.path.join(f"{tier_hit_path(setup)}", "{datatype}","{period}", "{run}", "{experiment}-{period}-{run}-{datatype}-{timestamp}-tier_hit.lh5")

def get_pattern_tier_evt(setup):
    return os.path.join(f"{tier_evt_path(setup)}", "{datatype}","{period}", "{run}", "{experiment}-{period}-{run}-{datatype}-{timestamp}-tier_evt.lh5")

def get_pattern_tier(setup,tier):
    if tier =="daq":
        return get_pattern_tier_daq(setup)
    elif tier =="raw":
        return get_pattern_tier_raw(setup)
    elif tier =="tcm":
        return get_pattern_tier_tcm(setup)
    elif tier =="dsp":
        return get_pattern_tier_dsp(setup)
    elif tier =="hit":
        return get_pattern_tier_hit(setup)
    elif tier =="evt":
        return get_pattern_tier_evt(setup)
    else:
        raise Exception("invalid tier")

def get_pattern_par_raw(setup, name=None):
    if name is not None:
        return os.path.join(f"{par_raw_path(setup)}",  "cal", "{period}", "{run}", "{experiment}-{period}-{run}-cal-{timestamp}-par_raw_"+name+".json")
    else:
        return os.path.join(f"{par_raw_path(setup)}", "cal", "{period}", "{run}", "{experiment}-{period}-{run}-cal-{timestamp}-par_raw.json")

def get_pattern_par_tcm(setup, name=None):
    if name is not None:
        return os.path.join(f"{par_tcm_path(setup)}",  "cal", "{period}", "{run}", "{experiment}-{period}-{run}-cal-{timestamp}-par_tcm_"+name+".json")
    else:
        return os.path.join(f"{par_tcm_path(setup)}",  "cal", "{period}", "{run}", "{experiment}-{period}-{run}-cal-{timestamp}-par_tcm.json")

def get_pattern_par_dsp(setup, name=None):
    if name is None:
        return os.path.join(f"{par_dsp_path(setup)}", "cal", "{period}", "{run}", "{experiment}-{period}-{run}-cal-{timestamp}-par_dsp.json")
    elif name == "energy_grid":
        return os.path.join(f"{par_dsp_path(setup)}", "cal",  "{period}",  "{run}", "{experiment}-{period}-{run}-cal-{timestamp}-par_dsp_energy_grid.pkl")
    else:
        return os.path.join(f"{par_dsp_path(setup)}",  "cal", "{period}", "{run}", "{experiment}-{period}-{run}-cal-{timestamp}-par_dsp_"+name+".json")

def get_pattern_par_hit(setup, name=None):
    if name is not None:
        return os.path.join(f"{par_hit_path(setup)}",  "cal", "{period}", "{run}", "{experiment}-{period}-{run}-cal-{timestamp}-par_hit_"+name+".json")
    else:
        return os.path.join(f"{par_hit_path(setup)}",  "cal", "{period}", "{run}", "{experiment}-{period}-{run}-cal-{timestamp}-par_hit.json")

def get_pattern_par_evt(setup, name=None):
    if name is not None:
        return os.path.join(f"{par_evt_path(setup)}", "cal", "{period}", "{run}", "{experiment}-{period}-{run}-cal-{timestamp}-par_evt_"+name+".json")
    else:
        return os.path.join(f"{par_evt_path(setup)}", "cal", "{period}", "{run}", "{experiment}-{period}-{run}-cal-{timestamp}-par_evt.json")


def get_pattern_pars(setup, tier, name = None):
    if tier =="raw":
        return get_pattern_par_raw(setup, name=name)
    elif tier =="tcm":
        return get_pattern_par_tcm(setup, name=name)
    elif tier =="dsp":
        return get_pattern_par_dsp(setup, name=name)
    elif tier =="hit":
        return get_pattern_par_hit(setup, name=name)
    elif tier =="evt":
        return get_pattern_par_evt(setup, name=name)
    else:
        raise Exception("invalid tier")

def get_pattern_pars_overwrite(setup, tier, name = None):
    if name is not None:
        return os.path.join(f"{par_overwrite_path(setup)}", tier,  "{datatype}", "{period}", "{run}", "{experiment}-{period}-{run}-{datatype}-{timestamp}-par_"+tier+"_"+name+"-overwrite.json")
    else:
        return os.path.join(f"{par_overwrite_path(setup)}", tier,  "{datatype}", "{period}", "{run}", "{experiment}-{period}-{run}-{datatype}-{timestamp}-par_"+tier+"-overwrite.json")

def get_pattern_pars_tmp_channel(setup, tier, name=None):
    if name is None:
        return os.path.join(f"{tmp_par_path(setup)}", "{experiment}-{period}-{run}-cal-{timestamp}-{channel}-par_"+tier+".json")
    elif name == "energy_grid":
        return os.path.join(f"{tmp_par_path(setup)}", "{experiment}-{period}-{run}-cal-{timestamp}-{channel}-par_dsp_energy_grid.pkl")
    else:
        return os.path.join(f"{tmp_par_path(setup)}", "{experiment}-{period}-{run}-cal-{timestamp}-{channel}-par_"+tier+"_"+name+".json")

def get_pattern_plts_tmp_channel(setup, tier, name=None):
    if name is None:
        return os.path.join(f"{tmp_plts_path(setup)}", "{experiment}-{period}-{run}-cal-{timestamp}-{channel}-plt_"+tier+".pkl")
    else:
        return os.path.join(f"{tmp_plts_path(setup)}", "{experiment}-{period}-{run}-cal-{timestamp}-{channel}-plt_"+tier+"_"+name+".pkl")

def get_pattern_plts(setup, tier):
    return os.path.join(f"{plts_path(setup)}", tier,"cal", "{period}", "{run}", "{experiment}-{period}-{run}-cal-{timestamp}-plt_"+tier+".pkl")

def get_energy_grids_pattern_combine(setup):
    return os.path.join(f"{tmp_par_path(setup)}", "dsp", "cal",  "{{period}}", "{{run}}" , "par_dsp_energy_grid", "{{channel}}", "{{experiment}}-{{period}}-{{run}}-cal-{{timestamp}}-{{channel}}-{peak}-par_dsp_energy_grid.pkl")

def get_pattern_log(setup, processing_step):
    return os.path.join(f"{log_path(setup)}",processing_step , "{experiment}-{period}-{run}-{datatype}-{timestamp}-"+processing_step+".log")

def get_pattern_log_channel(setup, processing_step):
    return os.path.join(f"{log_path(setup)}",processing_step , "{experiment}-{period}-{run}-cal-{timestamp}-{channel}-"+processing_step+".log")