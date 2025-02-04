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
from legenddataflow.execenv import execenv_smk_py_script


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
        f'{execenv_smk_py_script(config, "par_geds_dsp_tau")}'
        "--configs {configs} "
        "--log {log} "
        "--datatype {params.datatype} "
        "--timestamp {params.timestamp} "
        "--channel {params.channel} "
        "--plot_path {output.plots} "
        "--output_file {output.decay_const} "
        "--pulser_file {input.pulser} "
        "--raw_files {input.files}"


rule build_pars_evtsel_geds:
    input:
        files=os.path.join(
            filelist_path(config), "all-{experiment}-{period}-{run}-cal-raw.filelist"
        ),
        pulser_file=get_pattern_pars_tmp_channel(config, "tcm", "pulser_ids"),
        database=get_pattern_pars_tmp_channel(config, "dsp", "decay_constant"),
        raw_cal=get_blinding_curve_file,
    params:
        timestamp="{timestamp}",
        datatype="cal",
        channel="{channel}",
    output:
        peak_file=temp(get_pattern_pars_tmp_channel(config, "dsp", "peaks", "lh5")),
    log:
        get_pattern_log_channel(config, "par_dsp_event_selection", time),
    group:
        "par-dsp"
    resources:
        runtime=300,
        mem_swap=70,
    shell:
        f'{execenv_smk_py_script(config, "par_geds_dsp_evtsel")}'
        "--configs {configs} "
        "--log {log} "
        "--datatype {params.datatype} "
        "--timestamp {params.timestamp} "
        "--channel {params.channel} "
        "--peak_file {output.peak_file} "
        "--pulser_file {input.pulser_file} "
        "--decay_const {input.database} "
        "--raw_cal {input.raw_cal} "
        "--raw_filelist {input.files}"


# This rule builds the optimal energy filter parameters for the dsp using fft files
rule build_pars_dsp_nopt_geds:
    input:
        files=os.path.join(
            filelist_path(config), "all-{experiment}-{period}-{run}-fft-raw.filelist"
        ),
        database=get_pattern_pars_tmp_channel(config, "dsp", "decay_constant"),
        inplots=get_pattern_plts_tmp_channel(config, "dsp", "decay_constant"),
    params:
        timestamp="{timestamp}",
        datatype="cal",
        channel="{channel}",
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
        f'{execenv_smk_py_script(config, "par_geds_dsp_nopt")}'
        "--database {input.database} "
        "--configs {configs} "
        "--log {log} "
        "--datatype {params.datatype} "
        "--timestamp {params.timestamp} "
        "--channel {params.channel} "
        "--inplots {input.inplots} "
        "--plot_path {output.plots} "
        "--dsp_pars {output.dsp_pars_nopt} "
        "--raw_filelist {input.files}"


# This rule builds the dplms energy filter for the dsp using fft and cal files
rule build_pars_dsp_dplms_geds:
    input:
        fft_files=os.path.join(
            filelist_path(config), "all-{experiment}-{period}-{run}-fft-raw.filelist"
        ),
        peak_file=get_pattern_pars_tmp_channel(config, "dsp", "peaks", "lh5"),
        database=get_pattern_pars_tmp_channel(config, "dsp", "noise_optimization"),
        inplots=get_pattern_plts_tmp_channel(config, "dsp", "noise_optimization"),
    params:
        timestamp="{timestamp}",
        datatype="cal",
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
        f'{execenv_smk_py_script(config, "par_geds_dsp_dplms")}'
        "--fft_raw_filelist {input.fft_files} "
        "--peak_file {input.peak_file} "
        "--database {input.database} "
        "--inplots {input.inplots} "
        "--configs {configs} "
        "--log {log} "
        "--datatype {params.datatype} "
        "--timestamp {params.timestamp} "
        "--channel {params.channel} "
        "--dsp_pars {output.dsp_pars} "
        "--lh5_path {output.lh5_path} "
        "--plot_path {output.plots} "


# This rule builds the optimal energy filter parameters for the dsp using calibration dsp files
rule build_pars_dsp_eopt_geds:
    input:
        peak_file=get_pattern_pars_tmp_channel(config, "dsp", "peaks", "lh5"),
        decay_const=get_pattern_pars_tmp_channel(config, "dsp", "dplms"),
        inplots=get_pattern_plts_tmp_channel(config, "dsp", "dplms"),
    params:
        timestamp="{timestamp}",
        datatype="cal",
        channel="{channel}",
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
        f'{execenv_smk_py_script(config, "par_geds_dsp_eopt")}'
        "--log {log} "
        "--configs {configs} "
        "--datatype {params.datatype} "
        "--timestamp {params.timestamp} "
        "--channel {params.channel} "
        "--peak_file {input.peak_file} "
        "--inplots {input.inplots} "
        "--decay_const {input.decay_const} "
        "--plot_path {output.plots} "
        "--qbb_grid_path {output.qbb_grid} "
        "--final_dsp_pars {output.dsp_pars}"


rule build_svm_dsp_geds:
    input:
        hyperpars=lambda wildcards: get_input_par_file(
            wildcards, "dsp", "svm_hyperpars"
        ),
        train_data=lambda wildcards: str(
            get_input_par_file(wildcards, "dsp", "svm_hyperpars")
        ).replace("hyperpars.json", "train.lh5"),
    output:
        dsp_pars=get_pattern_pars(config, "dsp", "svm", "pkl"),
    log:
        str(get_pattern_log(config, "pars_dsp_svm", time)).replace("{datatype}", "cal"),
    group:
        "par-dsp-svm"
    resources:
        runtime=300,
    shell:
        f'{execenv_smk_py_script(config, "par_geds_dsp_svm_build")}'
        "--log {log} "
        "--train_data {input.train_data} "
        "--train_hyperpars {input.hyperpars} "
        "--output_file {output.dsp_pars}"


rule build_pars_dsp_svm_geds:
    input:
        dsp_pars=get_pattern_pars_tmp_channel(config, "dsp_eopt"),
        svm_file=get_pattern_pars(config, "dsp", "svm", "pkl"),
    output:
        dsp_pars=temp(get_pattern_pars_tmp_channel(config, "dsp")),
    log:
        get_pattern_log_channel(config, "pars_dsp_svm", time),
    group:
        "par-dsp"
    resources:
        runtime=300,
    shell:
        f'{execenv_smk_py_script(config, "par_geds_dsp_svm")}'
        "--log {log} "
        "--input_file {input.dsp_pars} "
        "--output_file {output.dsp_pars} "
        "--svm_file {input.svm_file}"
