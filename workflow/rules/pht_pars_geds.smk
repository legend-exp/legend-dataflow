"""
Snakemake rules for processing pht (partition hit) tier data. This is done in 4 steps:
- extraction of calibration curves(s) for each run for each channel from cal data
- extraction of psd calibration parameters and partition level energy fitting for each channel over whole partition from cal data
- combining of all channels into single pars files with associated plot and results files
- running build hit over all channels using par file
"""

from pathlib import Path
from legenddataflow.methods.create_pars_keylist import ParsKeyResolve, ParsCatalog
from legenddataflow.methods.paths import filelist_path, config_path, metadata_path
from legenddataflow.methods.patterns import (
    get_pattern_pars_tmp_channel,
    get_pattern_plts_tmp_channel,
    get_pattern_log_channel,
    get_pattern_plts,
    get_pattern_tier,
    get_pattern_pars_tmp,
    get_pattern_log,
    get_pattern_pars,
)
from legenddataflowscripts.workflow import execenv_pyexe, set_last_rule_name

pht_par_catalog = ParsKeyResolve.get_par_catalog(
    ["-*-*-*-cal"],
    get_pattern_tier(config, "raw", check_in_cycle=False),
    {"cal": ["par_pht"], "lar": ["par_pht"]},
)

intier = "psp"

qc_pht_rules = {}
partcal_pht_rules = {}
aoe_pht_rules = {}
lq_pht_rules = {}
for key, dataset in part.datasets.items():
    for partition in dataset.keys():
        input_files = part.get_filelists(partition, key, intier)
        pulser_files = [
            str(file).replace("par_pht", "par_tcm")
            for file in part.get_par_files(
                pht_par_catalog, partition, key, tier="pht", name="pulser_ids"
            )
        ]
        tstamp = part.get_timestamp(pht_par_catalog, partition, key, tier="pht")
        wildcard_constrain = part.get_wildcard_constraints(partition, key)

        # qc rule
        rule:
            input:
                cal_files=input_files,
                fft_files=part.get_filelists(partition, key, intier, datatype="fft"),
                pulser_files=pulser_files,
                overwrite_files=get_input_par_file(
                    config,
                    tier="pht",
                    timestamp=part.get_timestamp(
                        pht_par_catalog,
                        partition,
                        key,
                        tier="pht",
                    ),
                ),
            wildcard_constraints:
                channel=wildcard_constrain,
            params:
                datatype="cal",
                channel="{channel}" if key == "default" else key,
                timestamp=part.get_timestamp(
                    pht_par_catalog, partition, key, tier="pht"
                ),
                dsp_table_name=lambda wildcards: get_table_name(
                    channelmap_textdb,
                    config,
                    "cal",
                    part.get_timestamp(pht_par_catalog, partition, key, tier="pht"),
                    wildcards.channel,
                    "dsp",
                ),
                configs=config_path(config),
                meta=metadata_path(config),
            output:
                hit_pars=[
                    temp(file)
                    for file in part.get_par_files(
                        pht_par_catalog,
                        partition,
                        key,
                        tier="pht",
                        name="qc",
                    )
                ],
                plot_file=[
                    temp(file)
                    for file in part.get_plt_files(
                        pht_par_catalog,
                        partition,
                        key,
                        tier="pht",
                        name="qc",
                    )
                ],
            log:
                part.get_log_file(
                    pht_par_catalog,
                    partition,
                    key,
                    "pht",
                    time,
                    name="par_pht_qc",
                ),
            group:
                "par-pht"
            resources:
                mem_swap=len(part.get_filelists(partition, key, intier)) * 30,
                runtime=300,
            shell:
                execenv_pyexe(config, "par-geds-pht-qc") + "--log {log} "
                "--configs {params.configs} "
                "--metadata {params.meta} "
                "--datatype {params.datatype} "
                "--timestamp {params.timestamp} "
                "--channel {params.channel} "
                "--table-name {params.dsp_table_name} "
                "--save-path {output.hit_pars} "
                "--plot-path {output.plot_file} "
                "--overwrite-files {input.overwrite_files} "
                "--pulser-files {input.pulser_files} "
                "--fft-files {input.fft_files} "
                "--cal-files {input.cal_files}"

        set_last_rule_name(workflow, f"{key}-{partition}-build_pht_qc")

        if key in qc_pht_rules:
            qc_pht_rules[key].append(list(workflow.rules)[-1])
        else:
            qc_pht_rules[key] = [list(workflow.rules)[-1]]

            # partition level energy calibration
        rule:
            input:
                files=input_files,
                pulser_files=pulser_files,
                ecal_file=part.get_par_files(
                    pht_par_catalog,
                    partition,
                    key,
                    tier="pht",
                    name="energy_cal",
                ),
                eres_file=part.get_par_files(
                    pht_par_catalog,
                    partition,
                    key,
                    tier="pht",
                    name="energy_cal_objects",
                    extension="pkl",
                ),
                inplots=part.get_plt_files(
                    pht_par_catalog,
                    partition,
                    key,
                    tier="pht",
                    name="energy_cal",
                ),
            wildcard_constraints:
                channel=wildcard_constrain,
            params:
                datatype="cal",
                channel="{channel}" if key == "default" else key,
                timestamp=tstamp,
                dsp_table_name=lambda wildcards: get_table_name(
                    channelmap_textdb,
                    config,
                    "cal",
                    part.get_timestamp(pht_par_catalog, partition, key, tier="pht"),
                    wildcards.channel,
                    "dsp",
                ),
                configs=config_path(config),
                meta=metadata_path(config),
            output:
                hit_pars=[
                    temp(file)
                    for file in part.get_par_files(
                        pht_par_catalog,
                        partition,
                        key,
                        tier="pht",
                        name="partcal",
                    )
                ],
                partcal_results=[
                    temp(file)
                    for file in part.get_par_files(
                        pht_par_catalog,
                        partition,
                        key,
                        tier="pht",
                        name="partcal_objects",
                        extension="pkl",
                    )
                ],
                plot_file=[
                    temp(file)
                    for file in part.get_plt_files(
                        pht_par_catalog,
                        partition,
                        key,
                        tier="pht",
                        name="partcal",
                    )
                ],
            log:
                part.get_log_file(
                    pht_par_catalog,
                    partition,
                    key,
                    "pht",
                    time,
                    name="par_pht_partcal",
                ),
            group:
                "par-pht"
            resources:
                mem_swap=len(part.get_filelists(partition, key, intier)) * 15,
                runtime=300,
            shell:
                execenv_pyexe(config, "par-geds-pht-ecal-part") + "--log {log} "
                "--configs {params.configs} "
                "--metadata {params.meta} "
                "--datatype {params.datatype} "
                "--timestamp {params.timestamp} "
                "--inplots {input.inplots} "
                "--channel {params.channel} "
                "--table-name {params.dsp_table_name} "
                "--fit-results {output.partcal_results} "
                "--eres-file {input.eres_file} "
                "--hit-pars {output.hit_pars} "
                "--plot-file {output.plot_file} "
                "--ecal-file {input.ecal_file} "
                "--pulser-files {input.pulser_files} "
                "--input-files {input.files}"

        set_last_rule_name(
            workflow, f"{key}-{partition}-build_pht_energy_super_calibrations"
        )
        previous_rule = list(workflow.rules)[-1]
        if key in partcal_pht_rules:
            partcal_pht_rules[key].append(previous_rule)
        else:
            partcal_pht_rules[key] = [previous_rule]

            # aoe calibration
        rule:
            input:
                files=input_files,
                pulser_files=pulser_files,
                ecal_file=strip_channel_wildcard_constraint(
                    previous_rule.output.hit_pars
                ),
                eres_file=strip_channel_wildcard_constraint(
                    previous_rule.output.partcal_results
                ),
                inplots=strip_channel_wildcard_constraint(
                    previous_rule.output.plot_file
                ),
            wildcard_constraints:
                channel=wildcard_constrain,
            params:
                datatype="cal",
                channel="{channel}" if key == "default" else key,
                timestamp=tstamp,
                dsp_table_name=lambda wildcards: get_table_name(
                    channelmap_textdb,
                    config,
                    "cal",
                    part.get_timestamp(pht_par_catalog, partition, key, tier="pht"),
                    wildcards.channel,
                    "dsp",
                ),
                configs=config_path(config),
                meta=metadata_path(config),
            output:
                hit_pars=[
                    temp(file)
                    for file in part.get_par_files(
                        pht_par_catalog,
                        partition,
                        key,
                        tier="pht",
                        name="aoecal",
                    )
                ],
                aoe_results=[
                    temp(file)
                    for file in part.get_par_files(
                        pht_par_catalog,
                        partition,
                        key,
                        tier="pht",
                        name="aoecal_objects",
                        extension="pkl",
                    )
                ],
                plot_file=[
                    temp(file)
                    for file in part.get_plt_files(
                        pht_par_catalog,
                        partition,
                        key,
                        tier="pht",
                        name="aoecal",
                    )
                ],
            log:
                part.get_log_file(
                    pht_par_catalog,
                    partition,
                    key,
                    "pht",
                    time,
                    name="par_pht_aoe",
                ),
            group:
                "par-pht"
            resources:
                mem_swap=len(part.get_filelists(partition, key, intier)) * 15,
                runtime=300,
            shell:
                execenv_pyexe(config, "par-geds-pht-aoe") + "--log {log} "
                "--configs {params.configs} "
                "--metadata {params.meta} "
                "--datatype {params.datatype} "
                "--timestamp {params.timestamp} "
                "--inplots {input.inplots} "
                "--channel {params.channel} "
                "--table-name {params.dsp_table_name} "
                "--aoe-results {output.aoe_results} "
                "--eres-file {input.eres_file} "
                "--hit-pars {output.hit_pars} "
                "--plot-file {output.plot_file} "
                "--ecal-file {input.ecal_file} "
                "--pulser-files {input.pulser_files} "
                "--input-files {input.files}"

        set_last_rule_name(
            workflow, f"{key}-{partition}-build_pht_aoe_calibrations"
        )
        previous_rule = list(workflow.rules)[-1]
        if key in aoe_pht_rules:
            aoe_pht_rules[key].append(previous_rule)
        else:
            aoe_pht_rules[key] = [previous_rule]

            # lq calibration
        rule:
            input:
                files=input_files,
                pulser_files=pulser_files,
                ecal_file=strip_channel_wildcard_constraint(
                    previous_rule.output.hit_pars
                ),
                eres_file=strip_channel_wildcard_constraint(
                    previous_rule.output.aoe_results
                ),
                inplots=strip_channel_wildcard_constraint(
                    previous_rule.output.plot_file
                ),
            wildcard_constraints:
                channel=wildcard_constrain,
            params:
                datatype="cal",
                channel="{channel}" if key == "default" else key,
                timestamp=tstamp,
                dsp_table_name=lambda wildcards: get_table_name(
                    channelmap_textdb,
                    config,
                    "cal",
                    part.get_timestamp(pht_par_catalog, partition, key, tier="pht"),
                    wildcards.channel,
                    "dsp",
                ),
                configs=config_path(config),
                meta=metadata_path(config),
            output:
                hit_pars=[
                    temp(file)
                    for file in part.get_par_files(
                        pht_par_catalog,
                        partition,
                        key,
                        tier="pht",
                    )
                ],
                lq_results=[
                    temp(file)
                    for file in part.get_par_files(
                        pht_par_catalog,
                        partition,
                        key,
                        tier="pht",
                        name="objects",
                        extension="pkl",
                    )
                ],
                plot_file=[
                    temp(file)
                    for file in part.get_plt_files(
                        pht_par_catalog,
                        partition,
                        key,
                        tier="pht",
                    )
                ],
            log:
                part.get_log_file(
                    pht_par_catalog,
                    partition,
                    key,
                    "pht",
                    time,
                    name="par_pht_lq",
                ),
            group:
                "par-pht"
            resources:
                mem_swap=len(part.get_filelists(partition, key, intier)) * 15,
                runtime=300,
            shell:
                execenv_pyexe(config, "par-geds-pht-lq") + "--log {log} "
                "--configs {params.configs} "
                "--metadata {params.meta} "
                "--datatype {params.datatype} "
                "--timestamp {params.timestamp} "
                "--inplots {input.inplots} "
                "--channel {params.channel} "
                "--table-name {params.dsp_table_name} "
                "--lq-results {output.lq_results} "
                "--eres-file {input.eres_file} "
                "--hit-pars {output.hit_pars} "
                "--plot-file {output.plot_file} "
                "--ecal-file {input.ecal_file} "
                "--pulser-files {input.pulser_files} "
                "--input-files {input.files}"

        set_last_rule_name(workflow, f"{key}-{partition}-build_pht_lq_calibration")

        if key in lq_pht_rules:
            lq_pht_rules[key].append(list(workflow.rules)[-1])
        else:
            lq_pht_rules[key] = [list(workflow.rules)[-1]]


# Merged energy and a/e supercalibrations to reduce number of rules as they have same inputs/outputs
# This rule builds the a/e calibration using the calibration dsp files for the whole partition
rule build_pht_qc:
    input:
        cal_files=(
            Path(filelist_path(config))
            / f"all-{{experiment}}-{{period}}-{{run}}-cal-{intier}.filelist",
        ),
        fft_files=(
            Path(filelist_path(config))
            / f"all-{{experiment}}-{{period}}-{{run}}-cal-{intier}.filelist",
        ),
        pulser_files=get_pattern_pars_tmp_channel(config, "tcm", "pulser_ids"),
        overwrite_files=lambda wildcards: get_input_par_file(
            config, tier="pht", wildcards=wildcards
        ),
    params:
        datatype="cal",
        channel="{channel}",
        timestamp="{timestamp}",
        dsp_table_name=lambda wildcards: get_table_name(
            channelmap_textdb,
            config,
            "cal",
            wildcards.timestamp,
            wildcards.channel,
            "dsp",
        ),
        configs=config_path(config),
        meta=metadata_path(config),
    output:
        hit_pars=temp(get_pattern_pars_tmp_channel(config, "pht", "qc")),
        plot_file=temp(get_pattern_plts_tmp_channel(config, "pht", "qc")),
    log:
        get_pattern_log_channel(config, "par_pht_qc", time),
    group:
        "par-pht"
    resources:
        mem_swap=60,
        runtime=300,
    shell:
        execenv_pyexe(config, "par-geds-pht-qc") + "--log {log} "
        "--configs {params.configs} "
        "--metadata {params.meta} "
        "--datatype {params.datatype} "
        "--timestamp {params.timestamp} "
        "--channel {params.channel} "
        "--table-name {params.dsp_table_name} "
        "--save-path {output.hit_pars} "
        "--plot-path {output.plot_file} "
        "--overwrite-files {input.overwrite_files} "
        "--pulser-files {input.pulser_files} "
        "--fft-files {input.fft_files} "
        "--cal-files {input.cal_files}"


fallback_qc_rule = list(workflow.rules)[-1]

rule_order_list = []
ordered = OrderedDict(qc_pht_rules)
ordered.move_to_end("default")
for key, items in ordered.items():
    rule_order_list += [item.name for item in items]
rule_order_list.append(fallback_qc_rule.name)
workflow._ruleorder.add(*rule_order_list)  # [::-1]


# This rule builds the energy calibration using the calibration dsp files
rule build_per_energy_calibration:
    input:
        files=Path(filelist_path(config))
        / f"all-{{experiment}}-{{period}}-{{run}}-cal-{intier}.filelist",
        pulser=get_pattern_pars_tmp_channel(config, "tcm", "pulser_ids"),
        pht_dict=rules.build_pht_qc.output.hit_pars,
        inplots=rules.build_pht_qc.output.plot_file,
        ctc_dict=ancient(
            lambda wildcards: ParsCatalog.get_par_file(
                psp_par_catalog if intier == "psp" else dsp_par_catalog,
                config,
                wildcards.timestamp,
                intier,
            )
        ),
    params:
        timestamp="{timestamp}",
        datatype="cal",
        channel="{channel}",
        tier="pht",
        dsp_table_name=lambda wildcards: get_table_name(
            channelmap_textdb,
            config,
            "cal",
            wildcards.timestamp,
            wildcards.channel,
            "dsp",
        ),
        configs=config_path(config),
        meta=metadata_path(config),
    output:
        ecal_file=temp(get_pattern_pars_tmp_channel(config, "pht", "energy_cal")),
        results_file=temp(
            get_pattern_pars_tmp_channel(
                config, "pht", "energy_cal_objects", extension="pkl"
            )
        ),
        plot_file=temp(get_pattern_plts_tmp_channel(config, "pht", "energy_cal")),
    log:
        get_pattern_log_channel(config, "par_pht_energy_cal", time),
    group:
        "par-pht"
    resources:
        runtime=300,
    shell:
        execenv_pyexe(config, "par-geds-hit-ecal") + "--log {log} "
        "--datatype {params.datatype} "
        "--timestamp {params.timestamp} "
        "--channel {params.channel} "
        "--table-name {params.dsp_table_name} "
        "--configs {params.configs} "
        "--metadata {params.meta} "
        "--tier {params.tier} "
        "--plot-path {output.plot_file} "
        "--results-path {output.results_file} "
        "--save-path {output.ecal_file} "
        "--inplot-dict {input.inplots} "
        "--in-hit-dict {input.pht_dict} "
        "--ctc-dict {input.ctc_dict} "
        "--pulser-file {input.pulser} "
        "--files {input.files}"


# Merged energy and a/e supercalibrations to reduce number of rules as they have same inputs/outputs
# This rule builds the a/e calibration using the calibration dsp files for the whole partition
rule build_pht_energy_super_calibrations:
    input:
        files=Path(filelist_path(config))
        / f"all-{{experiment}}-{{period}}-{{run}}-cal-{intier}.filelist",
        pulser_files=get_pattern_pars_tmp_channel(config, "tcm", "pulser_ids"),
        ecal_file=rules.build_per_energy_calibration.output.ecal_file,
        eres_file=rules.build_per_energy_calibration.output.results_file,
        inplots=rules.build_per_energy_calibration.output.plot_file,
    params:
        datatype="cal",
        channel="{channel}",
        timestamp="{timestamp}",
        dsp_table_name=lambda wildcards: get_table_name(
            channelmap_textdb,
            config,
            "cal",
            wildcards.timestamp,
            wildcards.channel,
            "dsp",
        ),
        configs=config_path(config),
        meta=metadata_path(config),
    output:
        hit_pars=temp(get_pattern_pars_tmp_channel(config, "pht", "partcal")),
        partcal_results=temp(
            get_pattern_pars_tmp_channel(
                config, "pht", "partcal_objects", extension="pkl"
            )
        ),
        plot_file=temp(get_pattern_plts_tmp_channel(config, "pht", "partcal")),
    log:
        get_pattern_log_channel(config, "par_pht_partcal", time),
    group:
        "par-pht"
    resources:
        mem_swap=60,
        runtime=300,
    shell:
        execenv_pyexe(config, "par-geds-pht-ecal-part") + "--log {log} "
        "--configs {params.configs} "
        "--metadata {params.meta} "
        "--datatype {params.datatype} "
        "--timestamp {params.timestamp} "
        "--channel {params.channel} "
        "--table-name {params.dsp_table_name} "
        "--inplots {input.inplots} "
        "--fit-results {output.partcal_results} "
        "--eres-file {input.eres_file} "
        "--hit-pars {output.hit_pars} "
        "--plot-file {output.plot_file} "
        "--ecal-file {input.ecal_file} "
        "--pulser-files {input.pulser_files} "
        "--input-files {input.files}"


fallback_partcal_pht_rule = list(workflow.rules)[-1]

rule_order_list = []
ordered = OrderedDict(partcal_pht_rules)
ordered.move_to_end("default")
for key, items in ordered.items():
    rule_order_list += [item.name for item in items]
rule_order_list.append(fallback_partcal_pht_rule.name)
workflow._ruleorder.add(*rule_order_list)  # [::-1]


# Merged energy and a/e supercalibrations to reduce number of rules as they have same inputs/outputs
# This rule builds the a/e calibration using the calibration dsp files for the whole partition
rule build_pht_aoe_calibrations:
    input:
        files=Path(filelist_path(config))
        / f"all-{{experiment}}-{{period}}-{{run}}-cal-{intier}.filelist",
        pulser_files=get_pattern_pars_tmp_channel(config, "tcm", "pulser_ids"),
        ecal_file=rules.build_pht_energy_super_calibrations.output.hit_pars,
        eres_file=rules.build_pht_energy_super_calibrations.output.partcal_results,
        inplots=rules.build_pht_energy_super_calibrations.output.plot_file,
    params:
        datatype="cal",
        channel="{channel}",
        timestamp="{timestamp}",
        dsp_table_name=lambda wildcards: get_table_name(
            channelmap_textdb,
            config,
            "cal",
            wildcards.timestamp,
            wildcards.channel,
            "dsp",
        ),
        configs=config_path(config),
        meta=metadata_path(config),
    output:
        hit_pars=temp(get_pattern_pars_tmp_channel(config, "pht", "aoecal")),
        aoe_results=temp(
            get_pattern_pars_tmp_channel(
                config, "pht", "aoecal_objects", extension="pkl"
            )
        ),
        plot_file=temp(get_pattern_plts_tmp_channel(config, "pht", "aoecal")),
    log:
        get_pattern_log_channel(config, "par_pht_aoe_cal", time),
    group:
        "par-pht"
    resources:
        mem_swap=60,
        runtime=300,
    shell:
        execenv_pyexe(config, "par-geds-pht-aoe") + "--log {log} "
        "--configs {params.configs} "
        "--metadata {params.meta} "
        "--datatype {params.datatype} "
        "--timestamp {params.timestamp} "
        "--inplots {input.inplots} "
        "--channel {params.channel} "
        "--table-name {params.dsp_table_name} "
        "--aoe-results {output.aoe_results} "
        "--eres-file {input.eres_file} "
        "--hit-pars {output.hit_pars} "
        "--plot-file {output.plot_file} "
        "--ecal-file {input.ecal_file} "
        "--pulser-files {input.pulser_files} "
        "--input-files {input.files}"


fallback_aoecal_pht_rule = list(workflow.rules)[-1]

rule_order_list = []
ordered = OrderedDict(aoe_pht_rules)
ordered.move_to_end("default")
for key, items in ordered.items():
    rule_order_list += [item.name for item in items]
rule_order_list.append(fallback_aoecal_pht_rule.name)
workflow._ruleorder.add(*rule_order_list)  # [::-1]


# This rule builds the lq calibration using the calibration dsp files for the whole partition
rule build_pht_lq_calibration:
    input:
        files=Path(filelist_path(config))
        / f"all-{{experiment}}-{{period}}-{{run}}-cal-{intier}.filelist",
        pulser_files=get_pattern_pars_tmp_channel(config, "tcm", "pulser_ids"),
        ecal_file=rules.build_pht_aoe_calibrations.output.hit_pars,
        eres_file=rules.build_pht_aoe_calibrations.output.aoe_results,
        inplots=rules.build_pht_aoe_calibrations.output.plot_file,
    params:
        datatype="cal",
        channel="{channel}",
        timestamp="{timestamp}",
        dsp_table_name=lambda wildcards: get_table_name(
            channelmap_textdb,
            config,
            "cal",
            wildcards.timestamp,
            wildcards.channel,
            "dsp",
        ),
        configs=config_path(config),
        meta=metadata_path(config),
    output:
        hit_pars=temp(get_pattern_pars_tmp_channel(config, "pht")),
        lq_results=temp(
            get_pattern_pars_tmp_channel(config, "pht", "objects", extension="pkl")
        ),
        plot_file=temp(get_pattern_plts_tmp_channel(config, "pht")),
    log:
        get_pattern_log_channel(config, "par_pht_lq_cal", time),
    group:
        "par-pht"
    resources:
        mem_swap=60,
        runtime=300,
    shell:
        execenv_pyexe(config, "par-geds-pht-lq") + "--log {log} "
        "--configs {params.configs} "
        "--metadata {params.meta} "
        "--datatype {params.datatype} "
        "--timestamp {params.timestamp} "
        "--inplots {input.inplots} "
        "--channel {params.channel} "
        "--table-name {params.dsp_table_name} "
        "--lq-results {output.lq_results} "
        "--eres-file {input.eres_file} "
        "--hit-pars {output.hit_pars} "
        "--plot-file {output.plot_file} "
        "--ecal-file {input.ecal_file} "
        "--pulser-files {input.pulser_files} "
        "--input-files {input.files}"


fallback_lq_pht_rule = list(workflow.rules)[-1]

rule_order_list = []
ordered = OrderedDict(lq_pht_rules)
ordered.move_to_end("default")
for key, items in ordered.items():
    rule_order_list += [item.name for item in items]
rule_order_list.append(fallback_lq_pht_rule.name)
workflow._ruleorder.add(*rule_order_list)  # [::-1]
