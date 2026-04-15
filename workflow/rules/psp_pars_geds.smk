"""
Snakemake rules for processing psp (partition dsp) tier data.
- extraction of calibration curves(s) for each run for each channel from cal data
- extraction of psd calibration parameters and partition level energy fitting for each channel over whole partition from cal data
"""

from legenddataflow.methods import ParsKeyResolve
from legenddataflow.methods.patterns import (
    get_pattern_pars_tmp_channel,
    get_pattern_plts_tmp_channel,
    get_pattern_log_channel,
    get_pattern_log,
    get_pattern_pars,
    get_pattern_tier,
)
from legenddataflow.methods.paths import config_path
from legenddataflowscripts.workflow import execenv_pyexe, set_last_rule_name

# merge just needed rules here, eopt, nopt



psp_rules = {}
for key, dataset in part.datasets.items():
    for partition in dataset.keys():

        rule:
            input:
                dsp_pars=part.get_par_files(
                    dsp_par_catalog,
                    partition,
                    key,
                    tier="dsp",
                    name="eopt",
                ),
                dsp_objs=part.get_par_files(
                    dsp_par_catalog,
                    partition,
                    key,
                    tier="dsp",
                    name="objects",
                    extension="pkl",
                ),
                dsp_plots=part.get_plt_files(
                    dsp_par_catalog, partition, key, tier="dsp"
                ),
            output:
                psp_pars=temp(
                    part.get_par_files(
                        psp_par_catalog,
                        partition,
                        key,
                        tier="psp",
                        name="eopt",
                    )
                ),
                psp_objs=temp(
                    part.get_par_files(
                        psp_par_catalog,
                        partition,
                        key,
                        tier="psp",
                        name="objects",
                        extension="pkl",
                    )
                ),
                psp_plots=temp(
                    part.get_plt_files(
                        psp_par_catalog,
                        partition,
                        key,
                        tier="psp",
                    )
                ),
            log:
                part.get_log_file(
                    psp_par_catalog,
                    partition,
                    key,
                    "psp",
                    time,
                    name="par_psp",
                ),
            wildcard_constraints:
                channel=part.get_wildcard_constraints(partition, key),
            group:
                "par-psp"
            resources:
                runtime=300,
            params:
                datatype="cal",
                channel="{channel}" if key == "default" else key,
                timestamp=part.get_timestamp(
                    psp_par_catalog, partition, key, tier="psp"
                ),
                configs=ro(config_path(config)),
            shell:
                execenv_pyexe(config, "par-geds-psp-average") + "--log {log} "
                "--configs {params.configs} "
                "--datatype {params.datatype} "
                "--timestamp {params.timestamp} "
                "--channel {params.channel} "
                "--in-plots {input.dsp_plots} "
                "--out-plots {output.psp_plots} "
                "--in-obj {input.dsp_objs} "
                "--out-obj {output.psp_objs} "
                "--input {input.dsp_pars} "
                "--output {output.psp_pars} "

        set_last_rule_name(workflow, f"{key}-{partition}-build_par_psp")

        if key in psp_rules:
            psp_rules[key].append(list(workflow.rules)[-1])
        else:
            psp_rules[key] = [list(workflow.rules)[-1]]


        # This rule builds the dplms energy filter for the dsp using fft and cal files
        rule build_pars_psp_dplms_geds_fallback:
            input:
                fft_files=os.path.join(
                    filelist_path(config), "all-{experiment}-{period}-{run}-fft-raw.filelist"
                ),
                peak_file=rules.build_pars_evtsel_geds.output.peak_file,
                database=rules.build_pars_dsp_pz_geds.output.dsp_pars,
                inplots=rules.build_pars_dsp_pz_geds.output.plots,
            output:
                dsp_pars=temp(get_pattern_pars_tmp_channel(config, "dsp", "dplms")),
                lh5_path=temp(get_pattern_pars_tmp_channel(config, "dsp", extension="lh5")),
                plots=temp(get_pattern_plts_tmp_channel(config, "dsp", "dplms")),
            log:
                get_pattern_log_channel(config, "pars_psp_dplms", time),
            group:
                "par-psp"
            resources:
                runtime=300,
            params:
                config_file=lambda wildcards: get_config_files(
                    dataflow_configs_texdb,
                    wildcards.timestamp,
                    "cal",
                    wildcards.channel,
                    "pars_psp_dplms",
                    "dplms_pars",
                ),
                processing_chain=lambda wildcards: get_config_files(
                    dataflow_configs_texdb,
                    wildcards.timestamp,
                    "cal",
                    wildcards.channel,
                    "pars_psp_dplms",
                    "proc_chain",
                ),
                log_config=lambda wildcards: get_log_config(
                    dataflow_configs_texdb,
                    wildcards.timestamp,
                    "cal",
                    "pars_psp_dplms",
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
            shell:
                execenv_pyexe(config, "par-geds-psp-dplms") + "--peak-file {input.peak_file} "
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


# Merged energy and a/e supercalibrations to reduce number of rules as they have same inputs/outputs
# This rule builds the a/e calibration using the calibration dsp files for the whole partition
rule build_par_psp_fallback:
    input:
        dsp_pars=get_pattern_pars_tmp_channel(config, "dsp", "eopt"),
        dsp_objs=get_pattern_pars_tmp_channel(config, "dsp", "objects", extension="pkl"),
        dsp_plots=get_pattern_plts_tmp_channel(config, "dsp"),
    output:
        psp_pars=temp(get_pattern_pars_tmp_channel(config, "psp", "eopt")),
        psp_objs=temp(
            get_pattern_pars_tmp_channel(config, "psp", "objects", extension="pkl")
        ),
        psp_plots=temp(get_pattern_plts_tmp_channel(config, "psp")),
    log:
        get_pattern_log_channel(config, "pars_psp", time),
    group:
        "par-psp"
    resources:
        runtime=300,
    params:
        datatype="cal",
        channel="{channel}",
        timestamp="{timestamp}",
        configs=ro(config_path(config)),
    shell:
        execenv_pyexe(config, "par-geds-psp-average") + "--log {log} "
        "--configs {params.configs} "
        "--datatype {params.datatype} "
        "--timestamp {params.timestamp} "
        "--channel {params.channel} "
        "--in-plots {input.dsp_plots} "
        "--out-plots {output.psp_plots} "
        "--in-obj {input.dsp_objs} "
        "--out-obj {output.psp_objs} "
        "--input {input.dsp_pars} "
        "--output {output.psp_pars} "


fallback_psp_rule = list(workflow.rules)[-1]
rule_order_list = []
ordered = OrderedDict(psp_rules)
ordered.move_to_end("default")
for key, items in ordered.items():
    rule_order_list += [item.name for item in items]
rule_order_list.append(fallback_psp_rule.name)
workflow._ruleorder.add(*rule_order_list)  # [::-1]


# This rule builds the dplms energy filter for the dsp using fft and cal files
rule build_pars_psp_dplms_geds_fallback:
    input:
        fft_files=os.path.join(
            filelist_path(config), "all-{experiment}-{period}-{run}-fft-raw.filelist"
        ),
        peak_file=rules.build_pars_evtsel_geds.output.peak_file,
        database=rules.build_pars_dsp_pz_geds.output.dsp_pars,
        inplots=rules.build_pars_dsp_pz_geds.output.plots,
    output:
        dsp_pars=temp(get_pattern_pars_tmp_channel(config, "dsp", "dplms")),
        lh5_path=temp(get_pattern_pars_tmp_channel(config, "dsp", extension="lh5")),
        plots=temp(get_pattern_plts_tmp_channel(config, "dsp", "dplms")),
    log:
        get_pattern_log_channel(config, "pars_psp_dplms", time),
    group:
        "par-psp"
    resources:
        runtime=300,
    params:
        config_file=lambda wildcards: get_config_files(
            dataflow_configs_texdb,
            wildcards.timestamp,
            "cal",
            wildcards.channel,
            "pars_psp_dplms",
            "dplms_pars",
        ),
        processing_chain=lambda wildcards: get_config_files(
            dataflow_configs_texdb,
            wildcards.timestamp,
            "cal",
            wildcards.channel,
            "pars_psp_dplms",
            "proc_chain",
        ),
        log_config=lambda wildcards: get_log_config(
            dataflow_configs_texdb,
            wildcards.timestamp,
            "cal",
            "pars_psp_dplms",
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
    shell:
        execenv_pyexe(config, "par-geds-psp-dplms") + "--peak-file {input.peak_file} "
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


rule build_svm_psp:
    input:
        hyperpars=lambda wildcards: get_input_par_file(
            config=config, wildcards=wildcards, tier="psp", name="svm_hyperpars"
        ),
        train_data=lambda wildcards: str(
            get_input_par_file(
                config=config, wildcards=wildcards, tier="psp", name="svm_hyperpars"
            )
        ).replace("hyperpars.yaml", "train.lh5"),
    output:
        dsp_pars=get_pattern_pars(config, "psp", "svm", extension="pkl"),
    log:
        str(get_pattern_log(config, "pars_psp_svm", time)).replace("{datatype}", "cal"),
    group:
        "par-dsp-svm"
    resources:
        runtime=300,
    params:
        timestamp="{timestamp}",
        datatype="cal",
        configs=ro(config_path(config)),
    shell:
        execenv_pyexe(config, "par-geds-dsp-svm-build") + "--log {log} "
        "--train-data {input.train_data} "
        "--train-hyperpars {input.hyperpars} "
        "--output-file {output.dsp_pars} "
        "--timestamp {params.timestamp} "
        "--datatype {params.datatype} "
        "--configs {params.configs} "


rule build_pars_psp_svm:
    input:
        dsp_pars=get_pattern_pars_tmp_channel(config, "psp_eopt"),
        svm_model=get_pattern_pars(config, "psp", "svm", extension="pkl"),
    output:
        dsp_pars=temp(get_pattern_pars_tmp_channel(config, "psp")),
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
        "--svm-file {input.svm_model}"
