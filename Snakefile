import pathlib, os, json, sys
import scripts as ds
from scripts.util.patterns import *
from datetime import datetime

# Set with `snakemake --configfile=/path/to/your/config.json`
# configfile: "have/to/specify/path/to/your/config.json"

subst_vars_in_snakemake_config(workflow, config)

setup = config["setups"]["l200"]
configs = config_path(setup)
chan_maps = chan_map_path(setup)
swenv = runcmd(setup)

basedir = workflow.basedir

localrules: do_nothing, gen_filelist, autogen_output

rule do_nothing:
    input:

onstart:
    print("Starting workflow")
    shell(f'rm {pars_path(setup)}/validity.jsonl || true')
    ds.pars_key_resolve.write_par_catalog(setup,['-*-*-*-cal'], os.path.join(pars_path(setup),"dsp",'validity.jsonl'), 
                                        get_pattern_tier_raw(setup), {"cal":["par_dsp"], 'lar':['par_dsp']})
    ds.pars_key_resolve.write_par_catalog(setup,['-*-*-*-cal'], os.path.join(pars_path(setup),"hit",'validity.jsonl'), 
                                        get_pattern_tier_raw(setup), {"cal":["par_hit"], 'lar':['par_hit']})

onsuccess:
    print("Workflow finished, no error")
    shell("rm *.gen || true")
    shell(f'rm {filelist_path(setup)}/* || true')
    

#Placeholder, can email or maybe put message in slack
onerror:
    print("An error occurred :( ")

checkpoint gen_filelist:
    output:
        os.path.join(filelist_path(setup),"{label}-{tier}.{extension}list")
    params:
        setup = lambda wildcards: setup,
        search_pattern = lambda wildcards: get_pattern_tier_raw(setup),
        basedir = basedir,
        configs = configs,
        chan_maps = chan_maps
    script:
        "scripts/create_{wildcards.extension}list.py"


def read_filelist(wildcards):
    with checkpoints.gen_filelist.get(label=wildcards.label, tier=wildcards.tier, extension="file").output[0].open() as f:
        files = f.read().splitlines()
        return files 

rule gen_fileDB_config:
    output:
        "fdb_config.json"
    script:
        "scripts/gen_fiileDB_config.py"
       


# Create "{label}-{tier}.gen", based on "{label}.keylist" via
# "{label}-{tier}.filelist". Will implicitly trigger creation of all files
# in "{label}-{tier}.filelist".
# Example: "all[-{detector}[-{measurement}[-{run}[-{timestamp}]]]]-{tier}.gen":
rule autogen_output:
    input:
        filelist = read_filelist,
    output:
        gen_output = "{label}-{tier}.gen",
        summary_log = f"{log_path(setup)}/summary-"+"{label}-{tier}"+f"-{datetime.strftime(datetime.utcnow(), '%Y%m%dT%H%M%SZ')}.log",
        warning_log = f"{log_path(setup)}/warning-"+"{label}-{tier}"+f"-{datetime.strftime(datetime.utcnow(), '%Y%m%dT%H%M%SZ')}.log",
    params:
        log_path = tmp_log_path(setup), 
        tmp_par_path = os.path.join(tmp_par_path(setup), "*_db.json"),
        valid_keys_path = os.path.join(pars_path(setup), "valid_keys"),
        filedb_path = os.path.join(pars_path(setup), "filedb"),
        setup = lambda wildcards: setup,
        basedir = basedir
    script:
        "scripts/complete_run.py" 


#This rule builds the tcm files each raw file
rule build_tier_tcm:
    input:
        get_pattern_tier_raw(setup)
    output:
        get_pattern_tier_tcm(setup)
    log:
        get_pattern_log(setup, "tier_tcm")
    group: "tier-tcm"
    resources:
        runtime=300
    shell:
        "{swenv} python3 -B {basedir}/scripts/build_tcm.py --log {log} --configs {configs} {input} {output}"

def read_filelist_raw_cal_channel(wildcards):
    label = f"all-{wildcards.experiment}-{wildcards.period}-{wildcards.run}-cal"
    with checkpoints.gen_filelist.get(label=label, tier="raw", extension="file").output[0].open() as f:
        files = f.read().splitlines()
        return files


rule build_pars_dsp_tau:
    input:
        files = read_filelist_raw_cal_channel
    params:
        timestamp = "{timestamp}",
        datatype = "cal",
        channel = "{channel}"
    output:
        decay_const = temp(get_pattern_pars_tmp_channel(setup, "dsp", "decay_constant")),
        plots = temp(get_pattern_plts_tmp_channel(setup, "dsp","decay_constant")) 
    log:
        get_pattern_log_channel(setup, "par_dsp_decay_constant")
    group: "par-dsp"
    resources:
        runtime=300
    shell:
        "{swenv} python3 -B {basedir}/scripts/pars_dsp_tau.py --configs {configs} --log {log} --datatype {params.datatype} --timestamp {params.timestamp} --channel {params.channel} --plot_path {output.plots} --output_file {output.decay_const} {input.files} "


#This rule builds the optimal energy filter parameters for the dsp using calibration dsp files
rule build_pars_dsp_eopt:
    input:
        files = os.path.join(filelist_path(setup),"all-{experiment}-{period}-{run}-cal-raw.filelist"),
        decay_const = get_pattern_pars_tmp_channel(setup, "dsp", "decay_constant"),
        inplots = get_pattern_plts_tmp_channel(setup, "dsp","decay_constant")
    params:
        timestamp = "{timestamp}",
        datatype = "cal",
        channel = "{channel}"
    output:
        dsp_pars = temp(get_pattern_pars_tmp_channel(setup, "dsp")),
        qbb_grid = temp(get_pattern_pars_tmp_channel(setup, "dsp", "energy_grid")),
        plots = temp(get_pattern_plts_tmp_channel(setup, "dsp")) 
    log:
        get_pattern_log_channel(setup, "pars_dsp_eopt")
    group: "par-dsp"
    resources:
        runtime=300
    shell:
        "{swenv} python3 -B {basedir}/scripts/pars_dsp_eopt.py --log {log} --configs {configs}  --datatype {params.datatype} --timestamp {params.timestamp}  --channel {params.channel} --raw_filelist {input.files} --inplots {input.inplots} --decay_const {input.decay_const} --plot_path {output.plots} --qbb_grid_path {output.qbb_grid} --final_dsp_pars {output.dsp_pars}" # {input.peak_files}


def read_filelist_pars_dsp_cal_channel(wildcards):
    """
    This function will read the filelist of the channels and return a list of dsp files one for each channel
    """

    label=f"all-{wildcards.experiment}-{wildcards.period}-{wildcards.run}-cal-{wildcards.timestamp}-channels"
    with checkpoints.gen_filelist.get(label=label, tier="dsp", extension="chan").output[0].open() as f:
        files = f.read().splitlines()
        return files 

def read_filelist_plts_dsp_cal_channel(wildcards):
    """
    This function will read the filelist of the channels and return a list of dsp files one for each channel
    """

    label=f"all-{wildcards.experiment}-{wildcards.period}-{wildcards.run}-cal-{wildcards.timestamp}-channels"
    with checkpoints.gen_filelist.get(label=label, tier="dsp", extension="chan").output[0].open() as f:
        files = f.read().splitlines()
        files = [file.replace("par", "plt").replace("json", "pkl") for file in files]
        return files 

def read_filelist_pars_dsp_cal_channel_results(wildcards):
    """
    This function will read the filelist of the channels and return a list of dsp files one for each channel
    """
    label=f"all-{wildcards.experiment}-{wildcards.period}-{wildcards.run}-cal-{wildcards.timestamp}-channels"
    with checkpoints.gen_filelist.get(label=label, tier="dsp_energy_grid", extension="chan").output[0].open() as f:
        files = f.read().splitlines()
        return files 


rule build_pars_dsp:
    input:
        read_filelist_pars_dsp_cal_channel,
        read_filelist_pars_dsp_cal_channel_results,
        read_filelist_plts_dsp_cal_channel
    output:
        get_pattern_par_dsp(setup),
        get_pattern_par_dsp(setup, name="energy_grid"),
        get_pattern_plts(setup, "dsp")
    group: "merge-dsp"
    shell:
        "{swenv} python3 -B {basedir}/scripts/merge_channels.py --input {input} --output {output}"

def get_pars_dsp_file(wildcards):
    """
    This function will get the pars file for the run checking the pars_overwrite 
    """
    out = ds.pars_catalog.get_par_file(setup, wildcards.timestamp, "dsp")
    return out

rule build_dsp:
    input:
        raw_file = get_pattern_tier_raw(setup),
        tcm_file = get_pattern_tier_tcm(setup),
        pars_file = ancient(get_pars_dsp_file)
    params:
        timestamp = "{timestamp}",
        datatype = "{datatype}"
    output:
        tier_file = get_pattern_tier_dsp(setup),
        db_file = get_pattern_pars_tmp(setup, "dsp_db")
    log:
        get_pattern_log(setup, "tier_dsp")
    group: "tier-dsp"
    resources:
        runtime=300
    shell:
        "{swenv} python3 -B {basedir}/scripts/build_dsp.py --log {log} --configs {configs}  --datatype {params.datatype} --timestamp {params.timestamp} --input {input.raw_file} --output {output.tier_file} --db_file {output.db_file} --pars_file {input.pars_file} " 



def read_filelist_dsp_cal(wildcards):
    label = f"all-{wildcards.experiment}-{wildcards.period}-{wildcards.run}-cal"
    with checkpoints.gen_filelist.get(label=label, tier="dsp", extension="file").output[0].open() as f:
        files = f.read().splitlines()
        return files

#This rule builds the energy calibration using the calibration dsp files 
rule build_energy_calibration:
    input:
        files = read_filelist_dsp_cal,
        ctc_dict = ancient(get_pars_dsp_file)
    params:
        timestamp = "{timestamp}",
        datatype = "cal",
        channel = "{channel}"
    output:
        ecal_file = temp(get_pattern_pars_tmp_channel(setup, "hit", "energy_cal")),
        results_file = temp(get_pattern_pars_tmp_channel(setup, "hit", "energy_cal_results")),
        plot_file = temp(get_pattern_plts_tmp_channel(setup, "hit","energy_cal"))
    log:
        get_pattern_log_channel(setup, "pars_hit_energy_cal")
    group: "par-hit"
    resources:
        runtime=300
    shell:
        "{swenv} python3 -B {basedir}/scripts/pars_hit_ecal.py --log {log} --datatype {params.datatype} --timestamp {params.timestamp} --channel {params.channel} --configs {configs} --plot_path {output.plot_file} --save_path {output.ecal_file} --ctc_dict {input.ctc_dict} --results_path {output.results_file} --files {input.files}" # 


#This rule builds the a/e calibration using the calibration dsp files 
rule build_aoe_calibration:
    input:
        files = os.path.join(filelist_path(setup),"all-{experiment}-{period}-{run}-cal-dsp.filelist"),
        ecal_file = get_pattern_pars_tmp_channel(setup, "hit", "energy_cal"),
        eres_file = get_pattern_pars_tmp_channel(setup, "hit", "energy_cal_results"),
        inplots = get_pattern_plts_tmp_channel(setup, "hit","energy_cal")
    params:
        timestamp = "{timestamp}",
        datatype = "cal",
        channel = "{channel}"
    output:
        hit_pars = temp(get_pattern_pars_tmp_channel(setup, "hit")),
        aoe_results = temp(get_pattern_pars_tmp_channel(setup, "hit", "results")),
        plot_file = temp(get_pattern_plts_tmp_channel(setup, "hit"))
    log:
        get_pattern_log_channel(setup, "pars_hit_aoe_cal")
    group: "par-hit"
    resources:
        runtime=300
    shell:
        "{swenv} python3 -B {basedir}/scripts/pars_hit_aoe.py  --log {log} --configs {configs} --datatype {params.datatype} --timestamp {params.timestamp} --inplots {input.inplots} --channel {params.channel} --aoe_results {output.aoe_results} --hit_pars {output.hit_pars} --plot_file {output.plot_file} --eres_file {input.eres_file} --ecal_file {input.ecal_file} {input.files}"     


def read_filelist_pars_hit_cal_channel(wildcards):
    """
    This function will read the filelist of the channels and return a list of dsp files one for each channel
    """
    label=f"all-{wildcards.experiment}-{wildcards.period}-{wildcards.run}-cal-{wildcards.timestamp}-channels"
    with checkpoints.gen_filelist.get(label=label, tier="hit", extension="chan").output[0].open() as f:
        files = f.read().splitlines()
        return files 

def read_filelist_plts_hit_cal_channel(wildcards):
    """
    This function will read the filelist of the channels and return a list of dsp files one for each channel
    """

    label=f"all-{wildcards.experiment}-{wildcards.period}-{wildcards.run}-cal-{wildcards.timestamp}-channels"
    with checkpoints.gen_filelist.get(label=label, tier="hit", extension="chan").output[0].open() as f:
        files = f.read().splitlines()
        files = [file.replace("par", "plt").replace("json", "pkl") for file in files]
        return files 

def read_filelist_pars_hit_cal_channel_results(wildcards):
    """
    This function will read the filelist of the channels and return a list of dsp files one for each channel
    """
    label=f"all-{wildcards.experiment}-{wildcards.period}-{wildcards.run}-cal-{wildcards.timestamp}-channels"
    with checkpoints.gen_filelist.get(label=label, tier="hit_results", extension="chan").output[0].open() as f:
        files = f.read().splitlines()
        return files 




checkpoint build_pars_hit:
    input:
        read_filelist_pars_hit_cal_channel,
        read_filelist_pars_hit_cal_channel_results,
        read_filelist_plts_hit_cal_channel
    output:
        get_pattern_par_hit(setup),
        get_pattern_par_hit(setup, name="results"),
        get_pattern_plts(setup, "hit")
    group: "merge-hit"
    shell:
        "{swenv} python3 -B {basedir}/scripts/merge_channels.py --input {input} --output {output}"

def get_pars_hit_file(wildcards):
    """
    This function will get the pars file for the run checking the pars_overwrite 
    """
    return ds.pars_catalog.get_par_file(setup, wildcards.timestamp, "hit")

rule build_hit:
    input:
        dsp_file = get_pattern_tier_dsp(setup),
        pars_file = get_pars_hit_file
    output:
        tier_file = get_pattern_tier_hit(setup),
        db_file = get_pattern_pars_tmp(setup, "hit_db")
    params:
        timestamp = "{timestamp}",
        datatype = "{datatype}"
    log:
        get_pattern_log(setup, "tier_hit")
    group: "tier-hit"
    resources:
        runtime=300
    shell:
        "{swenv} python3 -B {basedir}/scripts/build_hit.py  --configs {configs} --log {log} --datatype {params.datatype} --timestamp {params.timestamp} --pars_file {input.pars_file} --output {output.tier_file} --input {input.dsp_file} --db_file {output.db_file}"
