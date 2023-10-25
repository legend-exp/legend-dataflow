"""
Snakefile for doing the higher stages of data processing (everything beyond build_raw)
This includes:
- building the tcm
- dsp parameter generation
- building dsp
- hit pars generation
- building hit
- building evt
- the same for partition level tiers
"""

import pathlib, os, json, sys
import scripts as ds
from scripts.util.pars_loading import pars_catalog
from scripts.util.patterns import *
from datetime import datetime
from collections import OrderedDict

# Set with `snakemake --configfile=/path/to/your/config.json`
# configfile: "have/to/specify/path/to/your/config.json"

subst_vars_in_snakemake_config(workflow, config)

setup = config["setups"]["l200"]
configs = config_path(setup)
chan_maps = chan_map_path(setup)
swenv = runcmd(setup)
part = ds.dataset_file(setup, os.path.join(configs, "partitions.json"))

basedir = workflow.basedir


localrules:
    gen_filelist,
    autogen_output,


ds.pars_key_resolve.write_par_catalog(
    ["-*-*-*-cal"],
    os.path.join(pars_path(setup), "pht", "validity.jsonl"),
    get_pattern_tier_raw(setup),
    {"cal": ["par_pht"], "lar": ["par_pht"]},
)


onstart:
    print("Starting workflow")
    shell(f"rm {pars_path(setup)}/dsp/validity.jsonl || true")
    shell(f"rm {pars_path(setup)}/hit/validity.jsonl || true")
    shell(f"rm {pars_path(setup)}/pht/validity.jsonl || true")
    shell(f"rm {pars_path(setup)}/raw/validity.jsonl || true")
    ds.pars_key_resolve.write_par_catalog(
        setup,
        ["-*-*-*-cal"],
        os.path.join(pars_path(setup), "raw", "validity.jsonl"),
        get_pattern_tier_raw(setup),
        {"cal": ["par_raw"]},
    )
    ds.pars_key_resolve.write_par_catalog(
        setup,
        ["-*-*-*-cal"],
        os.path.join(pars_path(setup), "dsp", "validity.jsonl"),
        get_pattern_tier_raw(setup),
        {"cal": ["par_dsp"], "lar": ["par_dsp"]},
    )
    ds.pars_key_resolve.write_par_catalog(
        setup,
        ["-*-*-*-cal"],
        os.path.join(pars_path(setup), "hit", "validity.jsonl"),
        get_pattern_tier_raw(setup),
        {"cal": ["par_hit"], "lar": ["par_hit"]},
    )
    ds.pars_key_resolve.write_par_catalog(
        setup,
        ["-*-*-*-cal"],
        os.path.join(pars_path(setup), "pht", "validity.jsonl"),
        get_pattern_tier_raw(setup),
        {"cal": ["par_pht"], "lar": ["par_pht"]},
    )


onsuccess:
    print("Workflow finished, no error")
    shell("rm *.gen || true")
    shell(f"rm {filelist_path(setup)}/* || true")


# Placeholder, can email or maybe put message in slack
onerror:
    print("An error occurred :( ")


checkpoint gen_filelist:
    output:
        os.path.join(filelist_path(setup), "{label}-{tier}.{extension}list"),
    params:
        setup=lambda wildcards: setup,
        search_pattern=lambda wildcards: get_pattern_tier_raw(setup),
        basedir=basedir,
        configs=configs,
        chan_maps=chan_maps,
        blinding=False,
    script:
        "scripts/create_{wildcards.extension}list.py"


def read_filelist(wildcards):
    with checkpoints.gen_filelist.get(
        label=wildcards.label, tier=wildcards.tier, extension="file"
    ).output[0].open() as f:
        files = f.read().splitlines()
        return files


rule gen_fileDB_config:
    output:
        "fdb_config.json",
    script:
        "scripts/gen_fiileDB_config.py"


# Create "{label}-{tier}.gen", based on "{label}.keylist" via
# "{label}-{tier}.filelist". Will implicitly trigger creation of all files
# in "{label}-{tier}.filelist".
# Example: "all[-{detector}[-{measurement}[-{run}[-{timestamp}]]]]-{tier}.gen":
rule autogen_output:
    input:
        filelist=read_filelist,
    output:
        gen_output="{label}-{tier}.gen",
        summary_log=f"{log_path(setup)}/summary-"
        + "{label}-{tier}"
        + f"-{datetime.strftime(datetime.utcnow(), '%Y%m%dT%H%M%SZ')}.log",
        warning_log=f"{log_path(setup)}/warning-"
        + "{label}-{tier}"
        + f"-{datetime.strftime(datetime.utcnow(), '%Y%m%dT%H%M%SZ')}.log",
    params:
        log_path=tmp_log_path(setup),
        tmp_par_path=os.path.join(tmp_par_path(setup), "*_db.json"),
        valid_keys_path=os.path.join(pars_path(setup), "valid_keys"),
        filedb_path=os.path.join(pars_path(setup), "filedb"),
        setup=lambda wildcards: setup,
        basedir=basedir,
    script:
        "scripts/complete_run.py"


# This rule builds the tcm files each raw file
rule build_tier_tcm:
    input:
        get_pattern_tier_raw(setup),
    output:
        get_pattern_tier_tcm(setup),
    log:
        get_pattern_log(setup, "tier_tcm"),
    group:
        "tier-tcm"
    resources:
        runtime=300,
    shell:
        "{swenv} python3 -B {basedir}/scripts/build_tcm.py --log {log} --configs {configs} {input} {output}"


def read_filelist_raw_cal_channel(wildcards):
    label = f"all-{wildcards.experiment}-{wildcards.period}-{wildcards.run}-cal"
    with checkpoints.gen_filelist.get(label=label, tier="raw", extension="file").output[
        0
    ].open() as f:
        files = f.read().splitlines()
        return files


rule build_pars_dsp_tau:
    input:
        files=read_filelist_raw_cal_channel,
    params:
        timestamp="{timestamp}",
        datatype="cal",
        channel="{channel}",
    output:
        decay_const=temp(get_pattern_pars_tmp_channel(setup, "dsp", "decay_constant")),
        plots=temp(get_pattern_plts_tmp_channel(setup, "dsp", "decay_constant")),
    log:
        get_pattern_log_channel(setup, "par_dsp_decay_constant"),
    group:
        "par-dsp"
    resources:
        runtime=300,
    shell:
        "{swenv} python3 -B {basedir}/scripts/pars_dsp_tau.py "
        "--configs {configs} "
        "--log {log} "
        "--datatype {params.datatype} "
        "--timestamp {params.timestamp} "
        "--channel {params.channel} "
        "--plot_path {output.plots} "
        "--output_file {output.decay_const} "
        "{input.files}"


# This rule builds the optimal energy filter parameters for the dsp using calibration dsp files
rule build_pars_dsp_eopt:
    input:
        files=os.path.join(
            filelist_path(setup), "all-{experiment}-{period}-{run}-cal-raw.filelist"
        ),
        decay_const=get_pattern_pars_tmp_channel(setup, "dsp", "decay_constant"),
        inplots=get_pattern_plts_tmp_channel(setup, "dsp", "decay_constant"),
    params:
        timestamp="{timestamp}",
        datatype="cal",
        channel="{channel}",
    output:
        dsp_pars=temp(get_pattern_pars_tmp_channel(setup, "dsp")),
        qbb_grid=temp(get_pattern_pars_tmp_channel(setup, "dsp", "energy_grid")),
        plots=temp(get_pattern_plts_tmp_channel(setup, "dsp")),
    log:
        get_pattern_log_channel(setup, "pars_dsp_eopt"),
    group:
        "par-dsp"
    resources:
        runtime=300,
    shell:
        "{swenv} python3 -B {basedir}/scripts/pars_dsp_eopt.py "
        "--log {log} "
        "--configs {configs} "
        "--datatype {params.datatype} "
        "--timestamp {params.timestamp} "
        "--channel {params.channel} "
        "--raw_filelist {input.files} "
        "--inplots {input.inplots} "
        "--decay_const {input.decay_const} "
        "--plot_path {output.plots} "
        "--qbb_grid_path {output.qbb_grid} "
        "--final_dsp_pars {output.dsp_pars}"


def read_filelist_pars_dsp_cal_channel(wildcards):
    """
    This function will read the filelist of the channels and return a list of dsp files one for each channel
    """

    label = f"all-{wildcards.experiment}-{wildcards.period}-{wildcards.run}-cal-{wildcards.timestamp}-channels"
    with checkpoints.gen_filelist.get(label=label, tier="dsp", extension="chan").output[
        0
    ].open() as f:
        files = f.read().splitlines()
        return files


def read_filelist_plts_dsp_cal_channel(wildcards):
    """
    This function will read the filelist of the channels and return a list of dsp files one for each channel
    """

    label = f"all-{wildcards.experiment}-{wildcards.period}-{wildcards.run}-cal-{wildcards.timestamp}-channels"
    with checkpoints.gen_filelist.get(label=label, tier="dsp", extension="chan").output[
        0
    ].open() as f:
        files = f.read().splitlines()
        files = [file.replace("par", "plt").replace("json", "pkl") for file in files]
        return files


def read_filelist_pars_dsp_cal_channel_results(wildcards):
    """
    This function will read the filelist of the channels and return a list of dsp files one for each channel
    """
    label = f"all-{wildcards.experiment}-{wildcards.period}-{wildcards.run}-cal-{wildcards.timestamp}-channels"
    with checkpoints.gen_filelist.get(
        label=label, tier="dsp_energy_grid", extension="chan"
    ).output[0].open() as f:
        files = f.read().splitlines()
        return files


rule build_pars_dsp:
    input:
        read_filelist_pars_dsp_cal_channel,
        read_filelist_pars_dsp_cal_channel_results,
        read_filelist_plts_dsp_cal_channel,
    output:
        get_pattern_par_dsp(setup),
        get_pattern_par_dsp(setup, name="energy_grid", extension="pkl"),
        get_pattern_plts(setup, "dsp"),
    group:
        "merge-dsp"
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
        raw_file=get_pattern_tier_raw(setup),
        tcm_file=get_pattern_tier_tcm(setup),
        pars_file=ancient(get_pars_dsp_file),
    params:
        timestamp="{timestamp}",
        datatype="{datatype}",
    output:
        tier_file=get_pattern_tier_dsp(setup),
        db_file=get_pattern_pars_tmp(setup, "dsp_db"),
    log:
        get_pattern_log(setup, "tier_dsp"),
    group:
        "tier-dsp"
    resources:
        runtime=300,
        mem_swap=30,
    shell:
        "{swenv} python3 -B {basedir}/scripts/build_dsp.py "
        "--log {log} "
        "--configs {configs} "
        "--datatype {params.datatype} "
        "--timestamp {params.timestamp} "
        "--input {input.raw_file} "
        "--output {output.tier_file} "
        "--db_file {output.db_file} "
        "--pars_file {input.pars_file}"


def read_filelist_dsp_cal(wildcards):
    label = f"all-{wildcards.experiment}-{wildcards.period}-{wildcards.run}-cal"
    with checkpoints.gen_filelist.get(label=label, tier="dsp", extension="file").output[
        0
    ].open() as f:
        files = f.read().splitlines()
        return files


def read_filelist_raw_cal(wildcards):
    label = f"all-{wildcards.experiment}-{wildcards.period}-{wildcards.run}-cal"
    with checkpoints.gen_filelist.get(label=label, tier="raw", extension="file").output[
        0
    ].open() as f:
        files = f.read().splitlines()
        return files


def get_blinding_curve(wildcards):
    # func to load in blinding curves
    par_files = ds.pars_catalog.get_calib_files(
        Path(par_overwrite_path(setup)) / "raw" / "validity.jsonl", wildcards.timestamp
    )
    if isinstance(par_files, str):
        return str(Path(par_overwrite_path(setup)) / "raw" / par_files)
    else:
        return [
            str(Path(par_overwrite_path(setup)) / "raw" / par_file)
            for par_file in par_files
        ]


rule build_blinding_check:
    """
    Runs a check on the daqenergy of the calibration run that the blinding curve given still applies,
    if so creates a file whose existence will be checked by the raw blinding before proceeding with blinding the phy data"""
    input:
        files=read_filelist_raw_cal,
        par_file=get_blinding_curve,
    params:
        timestamp="{timestamp}",
        datatype="cal",
        channel="{channel}",
    output:
        get_pattern_pars_tmp_channel(setup, "raw"),
    log:
        get_pattern_log_channel(setup, "pars_hit_blind_check"),
    group:
        "par-hit"
    resources:
        runtime=300,
    shell:
        "{swenv} python3 -B {basedir}/scripts/check_blinding.py "
        "--log {log} "
        "--datatype {params.datatype} "
        "--timestamp {params.timestamp} "
        "--channel {params.channel} "
        "--configs {configs} "
        "--output {output} "
        "--blind_curve {input.par_file} "
        "--files {input.files} "


# This rule builds the energy calibration using the calibration dsp files
rule build_energy_calibration:
    input:
        files=read_filelist_dsp_cal,
        ctc_dict=ancient(get_pars_dsp_file),
    params:
        timestamp="{timestamp}",
        datatype="cal",
        channel="{channel}",
    output:
        ecal_file=temp(get_pattern_pars_tmp_channel(setup, "hit", "energy_cal")),
        results_file=temp(
            get_pattern_pars_tmp_channel(
                setup, "hit", "energy_cal_results", extension="pkl"
            )
        ),
        plot_file=temp(get_pattern_plts_tmp_channel(setup, "hit", "energy_cal")),
    log:
        get_pattern_log_channel(setup, "pars_hit_energy_cal"),
    group:
        "par-hit"
    resources:
        runtime=300,
    shell:
        "{swenv} python3 -B {basedir}/scripts/pars_hit_ecal.py "
        "--log {log} "
        "--datatype {params.datatype} "
        "--timestamp {params.timestamp} "
        "--channel {params.channel} "
        "--configs {configs} "
        "--plot_path {output.plot_file} "
        "--results_path {output.results_file} "
        "--save_path {output.ecal_file} "
        "--ctc_dict {input.ctc_dict} "
        "--files {input.files}"


# This rule builds the a/e calibration using the calibration dsp files
rule build_aoe_calibration:
    input:
        files=os.path.join(
            filelist_path(setup), "all-{experiment}-{period}-{run}-cal-dsp.filelist"
        ),
        ecal_file=get_pattern_pars_tmp_channel(setup, "hit", "energy_cal"),
        eres_file=get_pattern_pars_tmp_channel(
            setup, "hit", "energy_cal_results", extension="pkl"
        ),
        inplots=get_pattern_plts_tmp_channel(setup, "hit", "energy_cal"),
    params:
        timestamp="{timestamp}",
        datatype="cal",
        channel="{channel}",
    output:
        hit_pars=temp(get_pattern_pars_tmp_channel(setup, "hit")),
        aoe_results=temp(
            get_pattern_pars_tmp_channel(setup, "hit", "results", extension="pkl")
        ),
        plot_file=temp(get_pattern_plts_tmp_channel(setup, "hit")),
    log:
        get_pattern_log_channel(setup, "pars_hit_aoe_cal"),
    group:
        "par-hit"
    resources:
        runtime=300,
    shell:
        "{swenv} python3 -B {basedir}/scripts/pars_hit_aoe.py "
        "--log {log} "
        "--configs {configs} "
        "--datatype {params.datatype} "
        "--timestamp {params.timestamp} "
        "--inplots {input.inplots} "
        "--channel {params.channel} "
        "--aoe_results {output.aoe_results} "
        "--eres_file {input.eres_file} "
        "--hit_pars {output.hit_pars} "
        "--plot_file {output.plot_file} "
        "--ecal_file {input.ecal_file} {input.files}"


def read_filelist_pars_hit_cal_channel(wildcards):
    """
    This function will read the filelist of the channels and return a list of dsp files one for each channel
    """
    label = f"all-{wildcards.experiment}-{wildcards.period}-{wildcards.run}-cal-{wildcards.timestamp}-channels"
    with checkpoints.gen_filelist.get(label=label, tier="hit", extension="chan").output[
        0
    ].open() as f:
        files = f.read().splitlines()
        return files


def read_filelist_plts_hit_cal_channel(wildcards):
    """
    This function will read the filelist of the channels and return a list of dsp files one for each channel
    """

    label = f"all-{wildcards.experiment}-{wildcards.period}-{wildcards.run}-cal-{wildcards.timestamp}-channels"
    with checkpoints.gen_filelist.get(label=label, tier="hit", extension="chan").output[
        0
    ].open() as f:
        files = f.read().splitlines()
        files = [file.replace("par", "plt").replace("json", "pkl") for file in files]
        return files


def read_filelist_pars_hit_cal_channel_results(wildcards):
    """
    This function will read the filelist of the channels and return a list of dsp files one for each channel
    """
    label = f"all-{wildcards.experiment}-{wildcards.period}-{wildcards.run}-cal-{wildcards.timestamp}-channels"
    with checkpoints.gen_filelist.get(
        label=label, tier="hit_results", extension="chan"
    ).output[0].open() as f:
        files = f.read().splitlines()
        files = [file.replace("json", "pkl") for file in files]
        return files


checkpoint build_pars_hit:
    input:
        read_filelist_pars_hit_cal_channel,
        read_filelist_pars_hit_cal_channel_results,
        read_filelist_plts_hit_cal_channel,
    output:
        get_pattern_par_hit(setup),
        get_pattern_par_hit(setup, name="results", extension="dir"),
        get_pattern_plts(setup, "hit"),
    group:
        "merge-hit"
    shell:
        "{swenv} python3 -B {basedir}/scripts/merge_channels.py --input {input} --output {output}"


def read_filelist_pars_raw_cal_channel(wildcards):
    """
    This function will read the filelist of the channels and return a list of dsp files one for each channel
    """
    label = f"all-{wildcards.experiment}-{wildcards.period}-{wildcards.run}-cal-{wildcards.timestamp}-channels"
    with checkpoints.gen_filelist.get(label=label, tier="raw", extension="chan").output[
        0
    ].open() as f:
        files = f.read().splitlines()
        return files


checkpoint build_pars_raw:
    input:
        read_filelist_pars_raw_cal_channel,
    output:
        get_pattern_par_raw(setup),
    group:
        "merge-blinding"
    shell:
        """
        {swenv} python3 -B {basedir}/scripts/merge_channels.py \
        --input {input} \
        --output {output}
        """


def get_pars_hit_file(wildcards):
    """
    This function will get the pars file for the run checking the pars_overwrite
    """
    return ds.pars_catalog.get_par_file(setup, wildcards.timestamp, "hit")


def get_pars_raw_file(wildcards):
    """
    This function will get the pars file for the run checking the pars_overwrite
    """
    return ds.pars_catalog.get_par_file(setup, wildcards.timestamp, "raw")


rule build_hit:
    input:
        dsp_file=get_pattern_tier_dsp(setup),
        pars_file=get_pars_hit_file,
        blind_check_file=get_pars_raw_file,
    output:
        tier_file=get_pattern_tier_hit(setup),
        db_file=get_pattern_pars_tmp(setup, "hit_db"),
    params:
        timestamp="{timestamp}",
        datatype="{datatype}",
        tier="hit",
    log:
        get_pattern_log(setup, "tier_hit"),
    group:
        "tier-hit"
    resources:
        runtime=300,
    shell:
        "{swenv} python3 -B {basedir}/scripts/build_hit.py "
        "--configs {configs} "
        "--log {log} "
        "--tier {params.tier} "
        "--datatype {params.datatype} "
        "--timestamp {params.timestamp} "
        "--pars_file {input.pars_file} "
        "--output {output.tier_file} "
        "--input {input.dsp_file} "
        "--db_file {output.db_file}"


# This rule builds the energy calibration using the calibration dsp files
rule build_per_energy_calibration:
    input:
        files=read_filelist_dsp_cal,
        ctc_dict=ancient(get_pars_dsp_file),
    params:
        timestamp="{timestamp}",
        datatype="cal",
        channel="{channel}",
    output:
        ecal_file=temp(get_pattern_pars_tmp_channel(setup, "pht", "energy_cal")),
        results_file=temp(
            get_pattern_pars_tmp_channel(
                setup, "pht", "energy_cal_results", extension="pkl"
            )
        ),
        plot_file=temp(get_pattern_plts_tmp_channel(setup, "pht", "energy_cal")),
    log:
        get_pattern_log_channel(setup, "pars_pht_energy_cal"),
    group:
        "par-pht"
    resources:
        runtime=300,
    shell:
        "{swenv} python3 -B {basedir}/scripts/pars_hit_ecal.py "
        "--log {log} "
        "--datatype {params.datatype} "
        "--timestamp {params.timestamp} "
        "--channel {params.channel} "
        "--configs {configs} "
        "--plot_path {output.plot_file} "
        "--results_path {output.results_file} "
        "--save_path {output.ecal_file} "
        "--ctc_dict {input.ctc_dict} "
        "--files {input.files}"


def read_filelist_pars_pht_cal_channel(wildcards):
    """
    This function will read the filelist of the channels and return a list of dsp files one for each channel
    """
    label = f"all-{wildcards.experiment}-{wildcards.period}-{wildcards.run}-cal-{wildcards.timestamp}-channels"
    with checkpoints.gen_filelist.get(label=label, tier="pht", extension="chan").output[
        0
    ].open() as f:
        files = f.read().splitlines()
        return files


def read_filelist_plts_pht_cal_channel(wildcards):
    """
    This function will read the filelist of the channels and return a list of dsp files one for each channel
    """

    label = f"all-{wildcards.experiment}-{wildcards.period}-{wildcards.run}-cal-{wildcards.timestamp}-channels"
    with checkpoints.gen_filelist.get(label=label, tier="pht", extension="chan").output[
        0
    ].open() as f:
        files = f.read().splitlines()
        files = [file.replace("par", "plt").replace("json", "pkl") for file in files]
        return files


def read_filelist_pars_pht_cal_channel_results(wildcards):
    """
    This function will read the filelist of the channels and return a list of dsp files one for each channel
    """
    label = f"all-{wildcards.experiment}-{wildcards.period}-{wildcards.run}-cal-{wildcards.timestamp}-channels"
    with checkpoints.gen_filelist.get(
        label=label, tier="pht_results", extension="chan"
    ).output[0].open() as f:
        files = f.read().splitlines()
        files = [file.replace("json", "pkl") for file in files]
        return files


checkpoint build_pars_pht:
    input:
        read_filelist_pars_pht_cal_channel,
        read_filelist_pars_pht_cal_channel_results,
        read_filelist_plts_pht_cal_channel,
    output:
        get_pattern_par_pht(setup),
        get_pattern_par_pht(setup, name="results", extension="dir"),
        get_pattern_plts(setup, "pht"),
    group:
        "merge-hit"
    shell:
        "{swenv} python3 -B {basedir}/scripts/merge_channels.py --input {input} --output {output}"


def get_pars_pht_file(wildcards):
    """
    This function will get the pars file for the run checking the pars_overwrite
    """
    return ds.pars_catalog.get_par_file(setup, wildcards.timestamp, "pht")


rule build_pht:
    input:
        dsp_file=get_pattern_tier_dsp(setup),
        #hit_file = get_pattern_tier_hit(setup),
        pars_file=get_pars_pht_file,
    output:
        tier_file=get_pattern_tier_pht(setup),
        db_file=get_pattern_pars_tmp(setup, "pht_db"),
    params:
        timestamp="{timestamp}",
        datatype="{datatype}",
        tier="pht",
    log:
        get_pattern_log(setup, "tier_pht"),
    group:
        "tier-pht"
    resources:
        runtime=300,
    shell:
        "{swenv} python3 -B {basedir}/scripts/build_hit.py "
        "--configs {configs} "
        "--log {log} "
        "--tier {params.tier} "
        "--datatype {params.datatype} "
        "--timestamp {params.timestamp} "
        "--pars_file {input.pars_file} "
        "--output {output.tier_file} "
        "--input {input.dsp_file} "
        "--db_file {output.db_file}"


# def fix_name(new_name):
#     """ sets the name of the most recently created rule to be new_name
#     """
#     list(workflow.rules)[-1].name = new_name
#     temp_rules = list(rules.__dict__.items())
#     temp_rules[-1] = (new_name, temp_rules[-1][1])
#     rules.__dict__ = dict(temp_rules)

part_pht_rules = {}
for key, dataset in part.datasets.items():
    for partition in dataset.keys():

        rule:
            input:
                files=part.get_filelists(partition, key, "dsp"),
                ecal_file=part.get_par_files(
                    f"{par_pht_path(setup)}/validity.jsonl",
                    partition,
                    key,
                    tier="pht",
                    name="energy_cal",
                ),
                eres_file=part.get_par_files(
                    f"{par_pht_path(setup)}/validity.jsonl",
                    partition,
                    key,
                    tier="pht",
                    name="energy_cal_results",
                    extension="pkl",
                ),
                inplots=part.get_plt_files(
                    f"{par_pht_path(setup)}/validity.jsonl",
                    partition,
                    key,
                    tier="pht",
                    name="energy_cal",
                ),
            params:
                datatype="cal",
                channel="{channel}",
                timestamp=part.get_timestamp(
                    f"{par_pht_path(setup)}/validity.jsonl", partition, key, tier="pht"
                ),
            output:
                hit_pars=[
                    temp(file)
                    for file in part.get_par_files(
                        f"{par_pht_path(setup)}/validity.jsonl",
                        partition,
                        key,
                        tier="pht",
                    )
                ],
                aoe_results=[
                    temp(file)
                    for file in part.get_par_files(
                        f"{par_pht_path(setup)}/validity.jsonl",
                        partition,
                        key,
                        tier="pht",
                        name="results",
                        extension="pkl",
                    )
                ],
                plot_file=[
                    temp(file)
                    for file in part.get_plt_files(
                        f"{par_pht_path(setup)}/validity.jsonl",
                        partition,
                        key,
                        tier="pht",
                    )
                ],
            log:
                part.get_log_file(
                    f"{par_pht_path(setup)}/validity.jsonl",
                    partition,
                    key,
                    "pht",
                    name="par_pht",
                ),
            group:
                "par-pht"
            resources:
                mem_swap=75,
                runtime=300,
            shell:
                "{swenv} python3 -B {basedir}/scripts/pars_pht.py "
                "--log {log} "
                "--configs {configs} "
                "--datatype {params.datatype} "
                "--timestamp {params.timestamp} "
                "--inplots {input.inplots} "
                "--channel {params.channel} "
                "--aoe_results {output.aoe_results} "
                "--eres_file {input.eres_file} "
                "--hit_pars {output.hit_pars} "
                "--plot_file {output.plot_file} "
                "--ecal_file {input.ecal_file} "
                "--input_files {input.files}"

        # fix_name(f"{key}-{partition}")

        if key in part_pht_rules:
            part_pht_rules[key].append(list(workflow.rules)[-1])
        else:
            part_pht_rules[key] = [list(workflow.rules)[-1]]


# Merged energy and a/e supercalibrations to reduce number of rules as they have same inputs/outputs
# This rule builds the a/e calibration using the calibration dsp files for the whole partition
rule build_pht_super_calibrations:
    input:
        files=os.path.join(
            filelist_path(setup), "all-{experiment}-{period}-{run}-cal-dsp.filelist"
        ),
        ecal_file=get_pattern_pars_tmp_channel(setup, "pht", "energy_cal"),
        eres_file=get_pattern_pars_tmp_channel(
            setup, "pht", "energy_cal_results", extension="pkl"
        ),
        inplots=get_pattern_plts_tmp_channel(setup, "pht", "energy_cal"),
    params:
        datatype="cal",
        channel="{channel}",
        timestamp="{timestamp}",
    output:
        hit_pars=temp(get_pattern_pars_tmp_channel(setup, "pht")),
        aoe_results=temp(
            get_pattern_pars_tmp_channel(setup, "pht", "results", extension="pkl")
        ),
        plot_file=temp(get_pattern_plts_tmp_channel(setup, "pht")),
    log:
        get_pattern_log_channel(setup, "pars_pht_aoe_cal"),
    group:
        "par-pht"
    resources:
        mem_swap=60,
        runtime=300,
    shell:
        "{swenv} python3 -B {basedir}/scripts/pars_pht.py "
        "--log {log} "
        "--configs {configs} "
        "--datatype {params.datatype} "
        "--timestamp {params.timestamp} "
        "--inplots {input.inplots} "
        "--channel {params.channel} "
        "--aoe_results {output.aoe_results} "
        "--eres_file {input.eres_file} "
        "--hit_pars {output.hit_pars} "
        "--plot_file {output.plot_file} "
        "--ecal_file {input.ecal_file} "
        "--input_files {input.files}"


fallback_pht_rule = list(workflow.rules)[-1]

rule_order_list = []
ordered = OrderedDict(part_pht_rules)
ordered.move_to_end("default")
for key, items in ordered.items():
    rule_order_list += [item.name for item in items]
rule_order_list.append(fallback_pht_rule.name)
workflow._ruleorder.add(*rule_order_list)  # [::-1]
