import os
from .utils import *

def key_pattern():
    return "{experiment}-{period}-{run}-{datatype}-{timestamp}"

def processing_pattern():
    return key_pattern()+'-{processing_step}'

def par_pattern():
    return "{experiment}-{period}-{run}-{datatype}-{timestamp}-pars_{name}"

def par_overwrite_pattern():
    return "{experiment}-{period}-{run}-{datatype}-{timestamp}-pars_{name}-overwrite"

def get_pattern_evts(setup, tier):
    if tier == "daq":
        return os.path.join(f"{inputdata_path(setup)}","daq", "{datatype}", "{period}", "{run}", "{experiment}-{period}-{run}-{datatype}-{timestamp}.fcio")
    else:
        return os.path.join(f"{evts_path(setup)}", tier, "{datatype}","{period}", "{run}", "{experiment}-{period}-{run}-{datatype}-{timestamp}-tier_" + tier + ".lh5")

def get_pattern_pars(setup, tier, name = None):
    if name is not None:
        return os.path.join(f"{pars_path(setup)}", tier,  "cal", "{period}", "{run}", "{experiment}-{period}-{run}-cal-{timestamp}-pars_"+tier+"_"+name+".json")
    else:
        return os.path.join(f"{pars_path(setup)}", tier,  "cal", "{period}", "{run}", "{experiment}-{period}-{run}-cal-{timestamp}-pars_"+tier+".json")

def get_pattern_pars_overwrite(setup, tier, name = None):
    if name is not None:
        return os.path.join(f"{par_overwrite_path(setup)}", tier,  "cal", "{period}", "{run}", "{experiment}-{period}-{run}-cal-{timestamp}-pars_"+tier+"_"+name+"-overwrite.json")
    else:
        return os.path.join(f"{par_overwrite_path(setup)}", tier,  "cal", "{period}", "{run}", "{experiment}-{period}-{run}-cal-{timestamp}-pars_"+tier+"-overwrite.json")

def get_pattern_pars_tmp_channel(setup, tier, name=None):
    if name =="energy_grid":
        return os.path.join(f"{tmp_par_path(setup)}", tier, "cal",  "{period}", "{run}" , "pars_dsp_energy_grid", "{channel}", "{experiment}-{period}-{run}-{channel}-{peak}-pars_dsp_energy_grid.pkl")
    elif name == "energy_grid_at_qbb":
        return os.path.join(f"{tmp_par_path(setup)}", tier, "cal",  "{period}",  "{run}", "pars_dsp_energy_grid", "{channel}", "{experiment}-{period}-{run}-{channel}-qbb-pars_dsp_energy_grid.pkl")
    else:
        if name is None:
            return os.path.join(f"{tmp_par_path(setup)}", tier, "cal",  "{period}", "{run}" ,  "{channel}" , "pars_"+tier, "{experiment}-{period}-{run}-{channel}-pars_"+tier+".json")
        else:
            return os.path.join(f"{tmp_par_path(setup)}", tier, "cal",  "{period}", "{run}" ,  "{channel}" , "pars_"+tier+"_"+name, "{experiment}-{period}-{run}-{channel}-pars_"+tier+"_"+name+".json")

def get_pattern_plts_tmp_channel(setup, tier, name):
    return os.path.join(f"{plts_path(setup)}", "plts",tier,"cal", "{period}", "{run}","{experiment}-{period}-{run}-{channel}-plts_"+tier+"_"+name+".pdf")

def get_pattern_plts(setup, tier, name):
    return os.path.join(f"{plts_path(setup)}", tier,"cal", "{period}", "{run}", "{experiment}-{period}-{run}-plts_"+tier+"_"+name+".pdf")

def get_energy_grids_pattern_combine(setup):
    return os.path.join(f"{tmp_par_path(setup)}", "dsp", "cal",  "{{period}}", "{{run}}" , "energy_grid", "{{channel}}", "{{experiment}}-{{period}}-{{run}}-{{channel}}-{peak}-pars_dsp_energy_grid.pkl")

