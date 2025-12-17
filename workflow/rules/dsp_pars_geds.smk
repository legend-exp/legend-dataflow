"""
Snakemake rules for building dsp pars for HPGes, before running build_dsp()
- extraction of pole zero constant(s) for each channel from cal data
- extraction of energy filter parameters and charge trapping correction for each channel from cal data
"""

from legenddataflow.methods import ParsKeyResolve
from legenddataflow.methods.patterns import (
    get_pattern_pars_tmp_channel,
    get_pattern_plts_tmp_channel,
    get_pattern_log_channel,
    get_pattern_tier,
    get_pattern_log,
    get_pattern_pars,
)
from legenddataflow.methods.paths import (
    filelist_path,
    config_path,
)
from legenddataflowscripts.workflow import execenv_pyexe


rule build_pars_dsp_pz_geds:
    input:
        files=(
            Path(filelist_path(config))
            / "all-{experiment}-{period}-{run}-cal-raw.filelist"
        ),
        pz_files=(
            Path(filelist_path(config))
            / "all-{experiment}-{period}-{run}-pzc-raw.filelist"
        ),
        pulser=get_pattern_pars_tmp_channel(config, "tcm", "pulser_ids"),
    params:
        config_file=lambda wildcards: get_config_files(
            dataflow_configs_texdb,
            wildcards.timestamp,
            "cal",
            wildcards.channel,
            "pars_dsp_tau",
            "tau_config",
        ),
        processing_chain=lambda wildcards: get_config_files(
            dataflow_configs_texdb,
            wildcards.timestamp,
            "cal",
            wildcards.channel,
            "pars_dsp_tau",
            "processing_chain",
        ),
        log_config=lambda wildcards: get_log_config(
            dataflow_configs_texdb,
            wildcards.timestamp,
            "cal",
            "pars_dsp_tau",
        ),
        raw_table_name=lambda wildcards: get_table_name(
            channelmap_textdb,
            config,
            "cal",
            wildcards.timestamp,
            wildcards.channel,
            "raw",
        ),
        configs=config_path(config),
    output:
        decay_const=temp(get_pattern_pars_tmp_channel(config, "dsp", "decay_constant")),
        plots=temp(get_pattern_plts_tmp_channel(config, "dsp", "decay_constant")),
    log:
        get_pattern_log_channel(config, "par_dsp_decay_constant", time),
    group:
        "par-dsp"
    resources:
        runtime=300,
    shell:
        execenv_pyexe(config, "par-geds-dsp-pz") + "--log {log} "
        "--log-config {params.log_config} "
        "--config-file {params.config_file} "
        "--processing-chain {params.processing_chain} "
        "--raw-table-name {params.raw_table_name} "
        "--plot-path {output.plots} "
        "--output-file {output.decay_const} "
        "--pulser-file {input.pulser} "
        "--raw-files {input.files} "
        "--pz-files {input.pz_files}"


rule build_pars_evtsel_geds:
    input:
        files=os.path.join(
            filelist_path(config), "all-{experiment}-{period}-{run}-cal-raw.filelist"
        ),
        pulser_file=get_pattern_pars_tmp_channel(config, "tcm", "pulser_ids"),
        database=rules.build_pars_dsp_pz_geds.output.decay_const,
        raw_cal_curve=get_blinding_curve_file,
    params:
        config_file=lambda wildcards: get_config_files(
            dataflow_configs_texdb,
            wildcards.timestamp,
            "cal",
            wildcards.channel,
            "pars_dsp_peak_selection",
            "peak_config",
        ),
        processing_chain=lambda wildcards: get_config_files(
            dataflow_configs_texdb,
            wildcards.timestamp,
            "cal",
            wildcards.channel,
            "pars_dsp_peak_selection",
            "processing_chain",
        ),
        log_config=lambda wildcards: get_log_config(
            dataflow_configs_texdb,
            wildcards.timestamp,
            "cal",
            "pars_dsp_peak_selection",
        ),
        raw_table_name=lambda wildcards: get_table_name(
            channelmap_textdb,
            config,
            "cal",
            wildcards.timestamp,
            wildcards.channel,
            "raw",
        ),
        configs=config_path(config),
        channel="{channel}",
    output:
        peak_file=temp(
            get_pattern_pars_tmp_channel(config, "dsp", "peaks", extension="lh5")
        ),
    log:
        get_pattern_log_channel(config, "par_dsp_event_selection", time),
    group:
        "par-dsp"
    resources:
        runtime=300,
        mem_swap=70,
    shell:
        execenv_pyexe(config, "par-geds-dsp-evtsel") + "--log {log} "
        "--log-config {params.log_config} "
        "--config-file {params.config_file} "
        "--processing-chain {params.processing_chain} "
        "--channel {params.channel} "
        "--raw-table-name {params.raw_table_name} "
        "--peak-file {output.peak_file} "
        "--pulser-file {input.pulser_file} "
        "--decay-const {input.database} "
        "--raw-cal-curve {input.raw_cal_curve} "
        "--raw-filelist {input.files}"


# This rule builds the optimal energy filter parameters for the dsp using fft files
rule build_pars_dsp_nopt_geds:
    input:
        files=os.path.join(
            filelist_path(config), "all-{experiment}-{period}-{run}-fft-raw.filelist"
        ),
        database=rules.build_pars_dsp_pz_geds.output.decay_const,
        inplots=rules.build_pars_dsp_pz_geds.output.plots,
    params:
        config_file=lambda wildcards: get_config_files(
            dataflow_configs_texdb,
            wildcards.timestamp,
            "cal",
            wildcards.channel,
            "pars_dsp_nopt",
            "optimiser_config",
        ),
        processing_chain=lambda wildcards: get_config_files(
            dataflow_configs_texdb,
            wildcards.timestamp,
            "cal",
            wildcards.channel,
            "pars_dsp_nopt",
            "processing_chain",
        ),
        log_config=lambda wildcards: get_log_config(
            dataflow_configs_texdb,
            wildcards.timestamp,
            "cal",
            "pars_dsp_nopt",
        ),
        raw_table_name=lambda wildcards: get_table_name(
            channelmap_textdb,
            config,
            "cal",
            wildcards.timestamp,
            wildcards.channel,
            "raw",
        ),
        configs=config_path(config),
    output:
        dsp_pars_nopt=temp(
            get_pattern_pars_tmp_channel(config, "dsp", "noise_optimization")
        ),
        plots=temp(get_pattern_plts_tmp_channel(config, "dsp", "noise_optimization")),
    log:
        get_pattern_log_channel(config, "par_dsp_noise_optimization", time),
    group:
        "par-dsp"
    resources:
        runtime=300,
    shell:
        execenv_pyexe(config, "par-geds-dsp-nopt") + "--database {input.database} "
        "--log {log} "
        "--log-config {params.log_config} "
        "--config-file {params.config_file} "
        "--processing-chain {params.processing_chain} "
        "--raw-table-name {params.raw_table_name} "
        "--inplots {input.inplots} "
        "--plot-path {output.plots} "
        "--dsp-pars {output.dsp_pars_nopt} "
        "--raw-filelist {input.files}"


# This rule builds the dplms energy filter for the dsp using fft and cal files
rule build_pars_dsp_dplms_geds:
    input:
        fft_files=os.path.join(
            filelist_path(config), "all-{experiment}-{period}-{run}-fft-raw.filelist"
        ),
        peak_file=rules.build_pars_evtsel_geds.output.peak_file,
        database=rules.build_pars_dsp_nopt_geds.output.dsp_pars_nopt,
        inplots=rules.build_pars_dsp_nopt_geds.output.plots,
    params:
        config_file=lambda wildcards: get_config_files(
            dataflow_configs_texdb,
            wildcards.timestamp,
            "cal",
            wildcards.channel,
            "pars_dsp_dplms",
            "dplms_pars",
        ),
        processing_chain=lambda wildcards: get_config_files(
            dataflow_configs_texdb,
            wildcards.timestamp,
            "cal",
            wildcards.channel,
            "pars_dsp_dplms",
            "proc_chain",
        ),
        log_config=lambda wildcards: get_log_config(
            dataflow_configs_texdb,
            wildcards.timestamp,
            "cal",
            "pars_dsp_dplms",
        ),
        raw_table_name=lambda wildcards: get_table_name(
            channelmap_textdb,
            config,
            "cal",
            wildcards.timestamp,
            wildcards.channel,
            "raw",
        ),
        configs=config_path(config),
        channel="{channel}",
    output:
        dsp_pars=temp(get_pattern_pars_tmp_channel(config, "dsp", "dplms")),
        lh5_path=temp(get_pattern_pars_tmp_channel(config, "dsp", extension="lh5")),
        plots=temp(get_pattern_plts_tmp_channel(config, "dsp", "dplms")),
    log:
        get_pattern_log_channel(config, "pars_dsp_dplms", time),
    group:
        "par-dsp"
    resources:
        runtime=300,
    shell:
        execenv_pyexe(config, "par-geds-dsp-dplms") + "--peak-file {input.peak_file} "
        "--fft-raw-filelist {input.fft_files} "
        "--database {input.database} "
        "--inplots {input.inplots} "
        "--log {log} "
        "--log-config {params.log_config} "
        "--config-file {params.config_file} "
        "--channel {params.channel} "
        "--processing-chain {params.processing_chain} "
        "--raw-table-name {params.raw_table_name} "
        "--dsp-pars {output.dsp_pars} "
        "--lh5-path {output.lh5_path} "
        "--plot-path {output.plots} "


# This rule builds the optimal energy filter parameters for the dsp using calibration dsp files
rule build_pars_dsp_eopt_geds:
    input:
        peak_file=rules.build_pars_evtsel_geds.output.peak_file,
        decay_const=rules.build_pars_dsp_dplms_geds.output.dsp_pars,
        inplots=rules.build_pars_dsp_dplms_geds.output.plots,
    params:
        config_file=lambda wildcards: get_config_files(
            dataflow_configs_texdb,
            wildcards.timestamp,
            "cal",
            wildcards.channel,
            "pars_dsp_eopt",
            "optimiser_config",
        ),
        processing_chain=lambda wildcards: get_config_files(
            dataflow_configs_texdb,
            wildcards.timestamp,
            "cal",
            wildcards.channel,
            "pars_dsp_eopt",
            "processing_chain",
        ),
        log_config=lambda wildcards: get_log_config(
            dataflow_configs_texdb,
            wildcards.timestamp,
            "cal",
            "pars_dsp_eopt",
        ),
        raw_table_name=lambda wildcards: get_table_name(
            channelmap_textdb,
            config,
            "cal",
            wildcards.timestamp,
            wildcards.channel,
            "raw",
        ),
    output:
        dsp_pars=temp(get_pattern_pars_tmp_channel(config, "dsp_eopt")),
        qbb_grid=temp(
            get_pattern_pars_tmp_channel(config, "dsp", "objects", extension="pkl")
        ),
        plots=temp(get_pattern_plts_tmp_channel(config, "dsp")),
    log:
        get_pattern_log_channel(config, "pars_dsp_eopt", time),
    group:
        "par-dsp"
    resources:
        runtime=300,
    shell:
        execenv_pyexe(config, "par-geds-dsp-eopt") + "--log {log} "
        "--log-config {params.log_config} "
        "--config-file {params.config_file} "
        "--processing-chain {params.processing_chain} "
        "--raw-table-name {params.raw_table_name} "
        "--peak-file {input.peak_file} "
        "--inplots {input.inplots} "
        "--decay-const {input.decay_const} "
        "--plot-path {output.plots} "
        "--qbb-grid-path {output.qbb_grid} "
        "--final-dsp-pars {output.dsp_pars}"


rule build_svm_dsp_geds:
    input:
        hyperpars=lambda wildcards: get_input_par_file(
            config=config,
            wildcards=wildcards,
            tier="dsp",
            name="svm_hyperpars",
            allow_none=True,
        ),
        train_data=lambda wildcards: (
            str(
                get_input_par_file(
                    config=config,
                    wildcards=wildcards,
                    tier="dsp",
                    name="svm_hyperpars",
                    allow_none=True,
                    overwrite=False,
                )
            ).replace("hyperpars.yaml", "train.lh5")
            if not isinstance(
                get_input_par_file(
                    config=config,
                    wildcards=wildcards,
                    tier="dsp",
                    name="svm_hyperpars",
                    allow_none=True,
                    overwrite=False,
                ),
                list,
            )
            else []
        ),
    params:
        timestamp="{timestamp}",
        datatype="cal",
        configs=config_path(config),
    output:
        dsp_pars=get_pattern_pars(config, "dsp", "svm", extension="pkl"),
    log:
        str(get_pattern_log(config, "pars_dsp_svm", time)).replace("{datatype}", "cal"),
    group:
        "par-dsp-svm"
    resources:
        runtime=300,
    shell:
        execenv_pyexe(config, "par-geds-dsp-svm-build") + "--log {log} "
        "--train-data {input.train_data} "
        "--train-hyperpars {input.hyperpars} "
        "--output-file {output.dsp_pars} "


rule build_pars_dsp_svm_geds:
    input:
        dsp_pars=rules.build_pars_dsp_eopt_geds.output.dsp_pars,
        svm_file=rules.build_svm_dsp_geds.output.dsp_pars,
    output:
        dsp_pars=temp(get_pattern_pars_tmp_channel(config, "dsp")),
    log:
        get_pattern_log_channel(config, "pars_dsp_svm", time),
    group:
        "par-dsp"
    resources:
        runtime=300,
    shell:
        execenv_pyexe(config, "par-geds-dsp-svm") + "--log {log} "
        "--input-file {input.dsp_pars} "
        "--output-file {output.dsp_pars} "
        "--svm-file {input.svm_file}"
