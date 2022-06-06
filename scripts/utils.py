import snakemake as smk
import re
import os
import shutil
import copy
import string

# For testing/debugging, use
# from scripts.utils import *
# import snakemake as smk
# setup = smk.load_configfile("config.json")["setups"]["l200"]

def inputdata_path(setup):
    return setup["paths"]["orig"]

def evts_path(setup):
    return setup["paths"]["gendata"]

def config_path(setup):
    return setup["paths"]["config"]

def pars_path(setup):
    return setup["paths"]["par"]

def tmp_par_path(setup):
    return setup["paths"]["tmp_par"]

def plts_path(setup):
    return setup["paths"]["plt"]

def par_overwrite_path(setup):
    return setup["paths"]["par_overwrite"]

def log_path(setup):
    return setup["paths"]["log"]


def key_pattern():
    return "{experiment}-{period}-{run}-{datatype}-{timestamp}"

def get_pattern_evts(setup, tier):
    if tier == "daq":
        return os.path.join(f"{inputdata_path(setup)}","daq", "{datatype}", "{period}", "{run}", "{experiment}-{period}-{run}-{datatype}-{timestamp}.fcio")
    else:
        return os.path.join(f"{evts_path(setup)}", tier, "{datatype}","{period}", "{run}", "{experiment}-{period}-{run}-{datatype}-{timestamp}-tier_" + tier + ".lh5")

def get_pattern_pars(setup, tier, name = None):
    if name is not None:
        return os.path.join(f"{pars_path(setup)}", tier,  "cal", "{period}", "{run}", "{experiment}-{period}-{run}-cal-pars_"+tier+"_"+name+".json")
    else:
        return os.path.join(f"{pars_path(setup)}", tier,  "cal", "{period}", "{run}", "{experiment}-{period}-{run}-cal-pars_"+tier+".json")

def get_pattern_pars_tmp_channel(setup, tier, name):
    if name =="energy_grid":
        return os.path.join(f"{tmp_par_path(setup)}", tier, "cal",  "{period}", "{run}" , "energy_grid", "{channel}", "{experiment}-{period}-{run}-{channel}-{peak}-pars_dsp_energy_grid.pkl")
    elif name == "energy_grid_at_qbb":
        return os.path.join(f"{tmp_par_path(setup)}", tier, "cal",  "{period}",  "{run}", "energy_grid_qbb", "{channel}", "{experiment}-{period}-{run}-{channel}-qbb-pars_dsp_energy_grid.pkl")
    else:
        return os.path.join(f"{tmp_par_path(setup)}", tier, "cal",  "{period}", "{run}" ,  "{channel}" , name, "{experiment}-{period}-{run}-{channel}-pars_"+tier+"_"+name+".json")

def get_pattern_plts_tmp_channel(setup, tier, name):
    return os.path.join(f"{plts_path(setup)}", "plts",tier,"cal", "{period}", "{run}","{experiment}-{period}-{run}-{channel}-plts_"+tier+"_"+name+".pdf")

def get_pattern_plts(setup, tier, name):
    return os.path.join(f"{plts_path(setup)}", tier,"cal", "{period}", "{run}", "{experiment}-{period}-{run}-plts_"+tier+"_"+name+".pdf")

def get_energy_grids_pattern_combine(setup):
    return os.path.join(f"{tmp_par_path(setup)}", "dsp", "cal",  "{{period}}", "{{run}}" , "energy_grid", "{{channel}}", "{{experiment}}-{{period}}-{{run}}-{{channel}}-{peak}-pars_dsp_energy_grid.pkl")

def runcmd(setup, envname):
    envcfg = setup["execenv"][envname]
    execcmd = envcfg["exec"]
    joined_cmd = execcmd if isinstance(execcmd, str) else " ".join(execcmd)
    envdefs = ""
    if "env" in envcfg:
        envdefs = " ".join([f"{key}={value}" for key, value in envcfg["env"].items()]) + " "
    return  envdefs + joined_cmd

def parse_keypart(keypart):
    keypart_rx = re.compile('(-(?P<experiment>[^-]+)(\\-(?P<period>[^-]+)(\\-(?P<run>[^-]+)(\\-(?P<datatype>[^-]+)(\\-(?P<timestamp>[^-]+))?)?)?)?)?$')
    d = keypart_rx.match(keypart).groupdict()
    for key in d:
        if d[key] is None:
            d[key] = "*"
    return d

def parse_channel_keypart(keypart):
    keypart_rx = re.compile('all(-(?P<experiment>[^-]+)(\\-(?P<period>[^-]+)(\\-(?P<run>[^-]+)(\\-(?P<datatype>[^-]+)(\\-(?P<channel>[^-]+))?)?)?)?)?$')
    d = keypart_rx.match(keypart).groupdict()
    for key in d:
        if d[key] is None:
            d[key] = "*"
    return d

def tier_files(setup, dataset_file, tier):
    key_pattern_rx = re.compile(smk.io.regex(key_pattern()))
    fn_pattern = get_pattern_evts(setup, tier)
    files = []
    with open(dataset_file) as f:
        for line in f:
            d = key_pattern_rx.match(line.strip()).groupdict()
            tier_filename = smk.io.expand(fn_pattern, experiment = d["experiment"], period = d["period"], run = d["run"], datatype  = d["datatype"], timestamp = d["timestamp"])[0]
            files.append(tier_filename)
    return files


def subst_vars_impl(x, var_values, ignore_missing = False):
    if isinstance(x, str):
        if '$' in x:
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

def subst_vars(props, var_values = {}, use_env = False, ignore_missing = False):
    combined_var_values = var_values
    if use_env:
        combined_var_values = {k: v for k, v in iter(os.environ.items())}
        combined_var_values.update(copy.copy(var_values))
    subst_vars_impl(props, combined_var_values, ignore_missing)


def subst_vars_in_snakemake_config(workflow, config):
    config_filename = workflow.overwrite_configfiles[0] # ToDo: Better way of handling this?
    subst_vars(
        config,
        var_values = {'_': os.path.dirname(config_filename)},
        use_env = True, ignore_missing = False
    )

def run_splitter(files):
    """
    Returns list containing lists of each run
    """
    
    runs = []
    run_files = []
    for file in files:
        base=os.path.basename(file)
        file_name = os.path.splitext(base)[0]
        parts = file_name.split('-')
        run_no = parts[3]
        if run_no not in runs:
            runs.append(run_no)
            run_files.append([])
        for i,run in enumerate(runs):
            if run == run_no:
                run_files[i].append(file) 
    return run_files
