"""
Snakemake rules for processing hit tier. This is done in 4 steps:
- extraction of calibration curves(s) for each channel from cal data
- extraction of psd calibration parameters for each channel from cal data
- combining of all channels into single pars files with associated plot and results files
- running build hit over all channels using par file
"""

from pathlib import Path
from legenddataflow.patterns import (
    get_pattern_pars_tmp_channel,
    get_pattern_plts_tmp_channel,
    get_pattern_log_channel,
    get_pattern_pars,
    get_pattern_plts,
    get_pattern_tier,
    get_pattern_pars_tmp,
    get_pattern_log,
    get_pattern_pars,
)
from legenddataflow.execenv import execenv_pyexe


# This rule builds the qc using the calibration dsp files and fft files
rule build_qc:
    input:
        files=os.path.join(
            filelist_path(config), "all-{experiment}-{period}-{run}-cal-dsp.filelist"
        ),
        fft_files=os.path.join(
            filelist_path(config), "all-{experiment}-{period}-{run}-fft-dsp.filelist"
        ),
        pulser=get_pattern_pars_tmp_channel(config, "tcm", "pulser_ids"),
        overwrite_files=lambda wildcards: get_overwrite_file("hit", wildcards),
    params:
        timestamp="{timestamp}",
        datatype="cal",
        channel="{channel}",
    output:
        qc_file=temp(get_pattern_pars_tmp_channel(config, "hit", "qc")),
        plot_file=temp(get_pattern_plts_tmp_channel(config, "hit", "qc")),
    log:
        get_pattern_log_channel(config, "pars_hit_qc", time),
    group:
        "par-hit"
    resources:
        runtime=300,
    shell:
        execenv_pyexe(config, "par-geds-hit-qc") + "--log {log} "
        "--datatype {params.datatype} "
        "--timestamp {params.timestamp} "
        "--channel {params.channel} "
        "--configs {configs} "
        "--metadata {meta} "
        "--plot-path {output.plot_file} "
        "--save-path {output.qc_file} "
        "--pulser-file {input.pulser} "
        "--cal-files {input.files} "
        "--fft-files {input.fft_files} "
        "--overwrite-files {input.overwrite_files} "


# This rule builds the energy calibration using the calibration dsp files
rule build_energy_calibration:
    input:
        files=os.path.join(
            filelist_path(config), "all-{experiment}-{period}-{run}-cal-dsp.filelist"
        ),
        pulser=get_pattern_pars_tmp_channel(config, "tcm", "pulser_ids"),
        ctc_dict=ancient(
            lambda wildcards: ParsCatalog.get_par_file(
                config, wildcards.timestamp, "dsp"
            )
        ),
        inplots=get_pattern_plts_tmp_channel(config, "hit", "qc"),
        in_hit_dict=get_pattern_pars_tmp_channel(config, "hit", "qc"),
    params:
        timestamp="{timestamp}",
        datatype="cal",
        channel="{channel}",
    output:
        ecal_file=temp(get_pattern_pars_tmp_channel(config, "hit", "energy_cal")),
        results_file=temp(
            get_pattern_pars_tmp_channel(
                config, "hit", "energy_cal_objects", extension="pkl"
            )
        ),
        plot_file=temp(get_pattern_plts_tmp_channel(config, "hit", "energy_cal")),
    log:
        get_pattern_log_channel(config, "pars_hit_energy_cal", time),
    group:
        "par-hit"
    resources:
        runtime=300,
    shell:
        execenv_pyexe(config, "par-geds-hit-ecal") + "--log {log} "
        "--datatype {params.datatype} "
        "--timestamp {params.timestamp} "
        "--channel {params.channel} "
        "--configs {configs} "
        "--metadata {meta} "
        "--plot-path {output.plot_file} "
        "--results-path {output.results_file} "
        "--save-path {output.ecal_file} "
        "--inplot-dict {input.inplots} "
        "--in-hit-dict {input.in_hit_dict} "
        "--ctc-dict {input.ctc_dict} "
        "--pulser-file {input.pulser} "
        "--files {input.files}"


# This rule builds the a/e calibration using the calibration dsp files
rule build_aoe_calibration:
    input:
        files=os.path.join(
            filelist_path(config), "all-{experiment}-{period}-{run}-cal-dsp.filelist"
        ),
        pulser=get_pattern_pars_tmp_channel(config, "tcm", "pulser_ids"),
        ecal_file=get_pattern_pars_tmp_channel(config, "hit", "energy_cal"),
        eres_file=get_pattern_pars_tmp_channel(
            config, "hit", "energy_cal_objects", extension="pkl"
        ),
        inplots=get_pattern_plts_tmp_channel(config, "hit", "energy_cal"),
    params:
        timestamp="{timestamp}",
        datatype="cal",
        channel="{channel}",
    output:
        hit_pars=temp(get_pattern_pars_tmp_channel(config, "hit", "aoe_cal")),
        aoe_results=temp(
            get_pattern_pars_tmp_channel(
                config, "hit", "aoe_cal_objects", extension="pkl"
            )
        ),
        plot_file=temp(get_pattern_plts_tmp_channel(config, "hit", "aoe_cal")),
    log:
        get_pattern_log_channel(config, "pars_hit_aoe_cal", time),
    group:
        "par-hit"
    resources:
        runtime=300,
    shell:
        execenv_pyexe(config, "par-geds-hit-aoe") + "--log {log} "
        "--configs {configs} "
        "--metadata {meta} "
        "--datatype {params.datatype} "
        "--timestamp {params.timestamp} "
        "--inplots {input.inplots} "
        "--channel {params.channel} "
        "--aoe-results {output.aoe_results} "
        "--eres-file {input.eres_file} "
        "--hit-pars {output.hit_pars} "
        "--plot-file {output.plot_file} "
        "--pulser-file {input.pulser} "
        "--ecal-file {input.ecal_file} "
        "{input.files}"


# This rule builds the lq calibration using the calibration dsp files
rule build_lq_calibration:
    input:
        files=os.path.join(
            filelist_path(config), "all-{experiment}-{period}-{run}-cal-dsp.filelist"
        ),
        pulser=get_pattern_pars_tmp_channel(config, "tcm", "pulser_ids"),
        ecal_file=get_pattern_pars_tmp_channel(config, "hit", "aoe_cal"),
        eres_file=get_pattern_pars_tmp_channel(
            config, "hit", "aoe_cal_objects", extension="pkl"
        ),
        inplots=get_pattern_plts_tmp_channel(config, "hit", "aoe_cal"),
    params:
        timestamp="{timestamp}",
        datatype="cal",
        channel="{channel}",
    output:
        hit_pars=temp(get_pattern_pars_tmp_channel(config, "hit")),
        lq_results=temp(
            get_pattern_pars_tmp_channel(config, "hit", "objects", extension="pkl")
        ),
        plot_file=temp(get_pattern_plts_tmp_channel(config, "hit")),
    log:
        get_pattern_log_channel(config, "pars_hit_lq_cal", time),
    group:
        "par-hit"
    resources:
        runtime=300,
    shell:
        execenv_pyexe(config, "par-geds-hit-lq") + "--log {log} "
        "--configs {configs} "
        "--metadata {meta} "
        "--datatype {params.datatype} "
        "--timestamp {params.timestamp} "
        "--inplots {input.inplots} "
        "--channel {params.channel} "
        "--lq-results {output.lq_results} "
        "--eres-file {input.eres_file} "
        "--hit-pars {output.hit_pars} "
        "--plot-file {output.plot_file} "
        "--pulser-file {input.pulser} "
        "--ecal-file {input.ecal_file} "
        "{input.files}"
