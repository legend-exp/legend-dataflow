"""
Snakemake rules for building dsp pars for HPGes, before running build_dsp()
- extraction of pole zero constant(s) for each channel from cal data
- extraction of energy filter parameters and charge trapping correction for each channel from cal data
"""

from legenddataflow.create_pars_keylist import ParsKeyResolve
from legenddataflow.patterns import (
    get_pattern_pars_tmp_channel,
    get_pattern_plts_tmp_channel,
    get_pattern_log_channel,
    get_pattern_tier,
    get_pattern_log,
    get_pattern_pars,
)
from legenddataflow.execenv import execenv_pyexe
from legenddataflow.scripts.table_name import get_table_name


rule build_pars_dsp_tau_geds:
    input:
        files=os.path.join(
            filelist_path(config), "all-{experiment}-{period}-{run}-cal-raw.filelist"
        ),
        pulser=get_pattern_pars_tmp_channel(config, "tcm", "pulser_ids"),
    params:
        timestamp="{timestamp}",
        datatype="cal",
        channel="{channel}",
        raw_table_name=lambda wildcards: get_table_name(
            channelmap_textdb,
            config,
            "cal",
            wildcards.timestamp,
            wildcards.channel,
            "raw",
        ),
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
        execenv_pyexe(config, "par-geds-dsp-tau") + "--configs {configs} "
        "--log {log} "
        "--datatype {params.datatype} "
        "--timestamp {params.timestamp} "
        "--channel {params.channel} "
        "--raw-table-name {params.raw_table_name} "
        "--plot-path {output.plots} "
        "--output-file {output.decay_const} "
        "--pulser-file {input.pulser} "
        "--raw-files {input.files}"


rule build_pars_evtsel_geds:
    input:
        files=os.path.join(
            filelist_path(config), "all-{experiment}-{period}-{run}-cal-raw.filelist"
        ),
        pulser_file=get_pattern_pars_tmp_channel(config, "tcm", "pulser_ids"),
        database=rules.build_pars_dsp_tau_geds.output.decay_const,
        raw_cal_curve=get_blinding_curve_file,
    params:
        timestamp="{timestamp}",
        datatype="cal",
        channel="{channel}",
        raw_table_name=lambda wildcards: get_table_name(
            channelmap_textdb,
            config,
            "cal",
            wildcards.timestamp,
            wildcards.channel,
            "raw",
        ),
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
        execenv_pyexe(config, "par-geds-dsp-evtsel") + "--configs {configs} "
        "--log {log} "
        "--datatype {params.datatype} "
        "--timestamp {params.timestamp} "
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
        database=rules.build_pars_dsp_tau_geds.output.decay_const,
        inplots=rules.build_pars_dsp_tau_geds.output.plots,
    params:
        timestamp="{timestamp}",
        datatype="cal",
        channel="{channel}",
        raw_table_name=lambda wildcards: get_table_name(
            channelmap_textdb,
            config,
            "cal",
            wildcards.timestamp,
            wildcards.channel,
            "raw",
        ),
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
        "--configs {configs} "
        "--log {log} "
        "--datatype {params.datatype} "
        "--timestamp {params.timestamp} "
        "--channel {params.channel} "
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
        timestamp="{timestamp}",
        datatype="cal",
        channel="{channel}",
        raw_table_name=lambda wildcards: get_table_name(
            channelmap_textdb,
            config,
            "cal",
            wildcards.timestamp,
            wildcards.channel,
            "raw",
        ),
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
        "--configs {configs} "
        "--log {log} "
        "--datatype {params.datatype} "
        "--timestamp {params.timestamp} "
        "--channel {params.channel} "
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
        timestamp="{timestamp}",
        datatype="cal",
        channel="{channel}",
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
        "--configs {configs} "
        "--datatype {params.datatype} "
        "--timestamp {params.timestamp} "
        "--channel {params.channel} "
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
            setup=config, wildcards=wildcards, tier="dsp", name="svm_hyperpars"
        ),
        train_data=lambda wildcards: str(
            get_input_par_file(
                setup=config, wildcards=wildcards, tier="dsp", name="svm_hyperpars"
            )
        ).replace("hyperpars.yaml", "train.lh5"),
    params:
        timestamp="{timestamp}",
        datatype="cal",
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
        "--timestamp {params.timestamp} "
        "--datatype {params.datatype} "
        "--configs {configs} "


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
