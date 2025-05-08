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
from legenddataflow.paths import config_path, metadata_path
from legenddataflow.execenv import execenv_pyexe


# This rule builds the qc using the calibration dsp files and fft files
rule par_hit_qc:
    input:
        files=os.path.join(
            filelist_path(config), "all-{experiment}-{period}-{run}-cal-dsp.filelist"
        ),
        fft_files=os.path.join(
            filelist_path(config), "all-{experiment}-{period}-{run}-fft-dsp.filelist"
        ),
        pulser=get_pattern_pars_tmp_channel(config, "tcm", "pulser_ids"),
        overwrite_files=lambda wildcards: get_input_par_file(
            config, tier="hit", wildcards=wildcards, allow_none=True
        ),
    params:
        timestamp="{timestamp}",
        datatype="cal",
        channel="{channel}",
        dsp_table_name=lambda wildcards: get_table_name(
            channelmap_textdb,
            config,
            "cal",
            wildcards.timestamp,
            wildcards.channel,
            "dsp",
        ),
        configs=ro(config_path(config)),
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
        "--configs {params.configs} "
        "--table-name {params.dsp_table_name} "
        "--plot-path {output.plot_file} "
        "--save-path {output.qc_file} "
        "--pulser-file {input.pulser} "
        "--cal-files {input.files} "
        "--fft-files {input.fft_files} "
        "--overwrite-files {input.overwrite_files} "


# This rule builds the energy calibration using the calibration dsp files
rule par_hit_energy_calibration:
    input:
        files=os.path.join(
            filelist_path(config), "all-{experiment}-{period}-{run}-cal-dsp.filelist"
        ),
        pulser=get_pattern_pars_tmp_channel(config, "tcm", "pulser_ids"),
        ctc_dict=ancient(
            lambda wildcards: ParsCatalog.get_par_file(
                dsp_par_catalog, config, wildcards.timestamp, "dsp"
            )
        ),
        inplots=rules.par_hit_qc.output.plot_file,
        in_hit_dict=rules.par_hit_qc.output.qc_file,
    params:
        timestamp="{timestamp}",
        datatype="cal",
        channel="{channel}",
        dsp_table_name=lambda wildcards: get_table_name(
            channelmap_textdb,
            config,
            "cal",
            wildcards.timestamp,
            wildcards.channel,
            "dsp",
        ),
        configs=ro(config_path(config)),
        meta=ro(metadata_path(config)),
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
        "--configs {params.configs} "
        "--metadata {params.meta} "
        "--table-name {params.dsp_table_name} "
        "--plot-path {output.plot_file} "
        "--results-path {output.results_file} "
        "--save-path {output.ecal_file} "
        "--inplot-dict {input.inplots} "
        "--in-hit-dict {input.in_hit_dict} "
        "--ctc-dict {input.ctc_dict} "
        "--pulser-file {input.pulser} "
        "--files {input.files}"


# This rule builds the a/e calibration using the calibration dsp files
rule par_hit_aoe_calibration:
    input:
        files=os.path.join(
            filelist_path(config), "all-{experiment}-{period}-{run}-cal-dsp.filelist"
        ),
        pulser=get_pattern_pars_tmp_channel(config, "tcm", "pulser_ids"),
        ecal_file=rules.par_hit_energy_calibration.output.ecal_file,
        eres_file=rules.par_hit_energy_calibration.output.results_file,
        inplots=rules.par_hit_energy_calibration.output.plot_file,
    params:
        timestamp="{timestamp}",
        datatype="cal",
        channel="{channel}",
        dsp_table_name=lambda wildcards: get_table_name(
            channelmap_textdb,
            config,
            "cal",
            wildcards.timestamp,
            wildcards.channel,
            "dsp",
        ),
        configs=ro(config_path(config)),
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
        "--configs {params.configs} "
        "--datatype {params.datatype} "
        "--timestamp {params.timestamp} "
        "--channel {params.channel} "
        "--table-name {params.dsp_table_name} "
        "--aoe-results {output.aoe_results} "
        "--hit-pars {output.hit_pars} "
        "--plot-file {output.plot_file} "
        "--inplots {input.inplots} "
        "--eres-file {input.eres_file} "
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
        ecal_file=rules.par_hit_aoe_calibration.output.hit_pars,
        eres_file=rules.par_hit_aoe_calibration.output.aoe_results,
        inplots=rules.par_hit_aoe_calibration.output.plot_file,
    params:
        timestamp="{timestamp}",
        datatype="cal",
        channel="{channel}",
        dsp_table_name=lambda wildcards: get_table_name(
            channelmap_textdb,
            config,
            "cal",
            wildcards.timestamp,
            wildcards.channel,
            "dsp",
        ),
        configs=ro(config_path(config)),
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
        "--configs {params.configs} "
        "--datatype {params.datatype} "
        "--timestamp {params.timestamp} "
        "--channel {params.channel} "
        "--table-name {params.dsp_table_name} "
        "--hit-pars {output.hit_pars} "
        "--plot-file {output.plot_file} "
        "--lq-results {output.lq_results} "
        "--pulser-file {input.pulser} "
        "--ecal-file {input.ecal_file} "
        "--inplots {input.inplots} "
        "--eres-file {input.eres_file} "
        "{input.files}"
