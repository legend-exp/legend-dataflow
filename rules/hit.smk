"""
Snakemake rules for processing hit tier. This is done in 4 steps:
- extraction of calibration curves(s) for each channel from cal data
- extraction of psd calibration parameters for each channel from cal data
- combining of all channels into single pars files with associated plot and results files
- running build hit over all channels using par file
"""

from scripts.util.pars_loading import pars_catalog
from scripts.util.patterns import (
    get_pattern_pars_tmp_channel,
    get_pattern_plts_tmp_channel,
    get_pattern_log_channel,
    get_pattern_par_hit,
    get_pattern_plts,
    get_pattern_tier_dsp,
    get_pattern_tier,
    get_pattern_pars_tmp,
    get_pattern_log,
    get_pattern_pars,
)


# This rule builds the energy calibration using the calibration dsp files
rule build_energy_calibration:
    input:
        files=lambda wildcards: read_filelist_cal(wildcards, "dsp"),
        tcm_filelist=os.path.join(
            filelist_path(setup), "all-{experiment}-{period}-{run}-cal-tcm.filelist"
        ),
        ctc_dict=ancient(
            lambda wildcards: pars_catalog.get_par_file(
                setup, wildcards.timestamp, "dsp"
            )
        ),
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
        "{swenv} python3 -B "
        f"{workflow.source_path('../scripts/pars_hit_ecal.py')} "
        "--log {log} "
        "--datatype {params.datatype} "
        "--timestamp {params.timestamp} "
        "--channel {params.channel} "
        "--configs {configs} "
        "--plot_path {output.plot_file} "
        "--results_path {output.results_file} "
        "--save_path {output.ecal_file} "
        "--ctc_dict {input.ctc_dict} "
        "--tcm_filelist {input.tcm_filelist} "
        "--files {input.files}"


# This rule builds the a/e calibration using the calibration dsp files
rule build_aoe_calibration:
    input:
        files=os.path.join(
            filelist_path(setup), "all-{experiment}-{period}-{run}-cal-dsp.filelist"
        ),
        tcm_filelist=os.path.join(
            filelist_path(setup), "all-{experiment}-{period}-{run}-cal-tcm.filelist"
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
        "{swenv} python3 -B "
        f"{workflow.source_path('../scripts/pars_hit_aoe.py')} "
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
        "--tcm_filelist {input.tcm_filelist} "
        "--ecal_file {input.ecal_file} "
        "{input.files}"


checkpoint build_pars_hit:
    input:
        lambda wildcards: read_filelist_pars_cal_channel(wildcards, "hit"),
        lambda wildcards: read_filelist_plts_cal_channel(wildcards, "hit"),
        lambda wildcards: read_filelist_pars_cal_channel(wildcards, "hit_results"),
    output:
        get_pattern_pars(setup, "hit", check_in_cycle=check_in_cycle),
        get_pattern_pars(
            setup,
            "hit",
            name="results",
            extension="dir",
            check_in_cycle=check_in_cycle,
        ),
        get_pattern_plts(setup, "hit"),
    group:
        "merge-hit"
    shell:
        "{swenv} python3 -B "
        f"{workflow.source_path('../scripts/merge_channels.py')} "
        "--input {input} "
        "--output {output} "


rule build_hit:
    input:
        dsp_file=get_pattern_tier_dsp(setup),
        pars_file=lambda wildcards: pars_catalog.get_par_file(
            setup, wildcards.timestamp, "hit"
        ),
    output:
        tier_file=get_pattern_tier(setup, "hit", check_in_cycle=check_in_cycle),
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
        "{swenv} python3 -B "
        f"{workflow.source_path('../scripts/build_hit.py')} "
        "--configs {configs} "
        "--log {log} "
        "--tier {params.tier} "
        "--datatype {params.datatype} "
        "--timestamp {params.timestamp} "
        "--pars_file {input.pars_file} "
        "--output {output.tier_file} "
        "--input {input.dsp_file} "
        "--db_file {output.db_file}"
