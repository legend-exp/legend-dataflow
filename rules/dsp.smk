"""
Snakemake rules for processing dsp tier. This is done in 4 steps:
- extraction of pole zero constant(s) for each channel from cal data
- extraction of energy filter parameters and charge trapping correction for each channel from cal data
- combining of all channels into single pars files with associated plot and results files
- running dsp over all channels using par file
"""
from scripts.util.pars_loading import pars_catalog
from scripts.util.patterns import (
    get_pattern_pars_tmp_channel,
    get_pattern_plts_tmp_channel,
    get_pattern_log_channel,
    get_pattern_par_dsp,
    get_pattern_plts,
    get_pattern_tier_raw,
    get_pattern_tier_tcm,
    get_pattern_tier_dsp,
    get_pattern_pars_tmp,
    get_pattern_log,
)


rule build_pars_dsp_tau:
    input:
        files=os.path.join(
            filelist_path(setup), "all-{experiment}-{period}-{run}-cal-raw.filelist"
        ),
        tcm_files=lambda wildcards: read_filelist_cal(wildcards, "tcm"),
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
        "{swenv} python3 -B "
        f"{workflow.source_path('../scripts/pars_dsp_tau.py')} "
        "--configs {configs} "
        "--log {log} "
        "--datatype {params.datatype} "
        "--timestamp {params.timestamp} "
        "--channel {params.channel} "
        "--plot_path {output.plots} "
        "--output_file {output.decay_const} "
        "--tcm_files {input.tcm_files} "
        "--raw_files {input.files}"


# This rule builds the optimal energy filter parameters for the dsp using calibration dsp files
rule build_pars_dsp_eopt:
    input:
        files=os.path.join(
            filelist_path(setup), "all-{experiment}-{period}-{run}-cal-raw.filelist"
        ),
        tcm_filelist=os.path.join(
            filelist_path(setup), "all-{experiment}-{period}-{run}-cal-tcm.filelist"
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
        "{swenv} python3 -B "
        f"{workflow.source_path('../scripts/pars_dsp_eopt.py')} "
        "--log {log} "
        "--configs {configs} "
        "--datatype {params.datatype} "
        "--timestamp {params.timestamp} "
        "--channel {params.channel} "
        "--raw_filelist {input.files} "
        "--tcm_filelist {input.tcm_filelist} "
        "--inplots {input.inplots} "
        "--decay_const {input.decay_const} "
        "--plot_path {output.plots} "
        "--qbb_grid_path {output.qbb_grid} "
        "--final_dsp_pars {output.dsp_pars}"


rule build_pars_dsp:
    input:
        lambda wildcards: read_filelist_pars_cal_channel(wildcards, "dsp"),
        lambda wildcards: read_filelist_plts_cal_channel(wildcards, "dsp"),
        lambda wildcards: read_filelist_pars_cal_channel(wildcards, "dsp_results"),
    output:
        get_pattern_par_dsp(setup),
        get_pattern_par_dsp(setup, name="results", extension="pkl"),
        get_pattern_plts(setup, "dsp"),
    group:
        "merge-dsp"
    shell:
        "{swenv} python3 -B "
        f"{workflow.source_path('../scripts/merge_channels.py')} "
        "--input {input} "
        "--output {output} "


# rule build_dsp:
#     input:
#         raw_file=get_pattern_tier_raw(setup),
#         tcm_file=get_pattern_tier_tcm(setup),
#         pars_file=ancient(
#             lambda wildcards: pars_catalog.get_par_file(
#                 setup, wildcards.timestamp, "dsp"
#             )
#         ),
#     params:
#         timestamp="{timestamp}",
#         datatype="{datatype}",
#     output:
#         tier_file=get_pattern_tier_dsp(setup),
#         db_file=get_pattern_pars_tmp(setup, "dsp_db"),
#     log:
#         get_pattern_log(setup, "tier_dsp"),
#     group:
#         "tier-dsp"
#     resources:
#         runtime=300,
#         mem_swap=30,
#     shell:
#         "{swenv} python3 -B "
#         f"{workflow.source_path('../scripts/build_dsp.py')} "
#         "--log {log} "
#         "--configs {configs} "
#         "--datatype {params.datatype} "
#         "--timestamp {params.timestamp} "
#         "--input {input.raw_file} "
#         "--output {output.tier_file} "
#         "--db_file {output.db_file} "
#         "--pars_file {input.pars_file}"
