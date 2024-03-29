"""
Snakemake rules for processing pht (partition hit) tier data. This is done in 4 steps:
- extraction of calibration curves(s) for each run for each channel from cal data
- extraction of psd calibration parameters and partition level energy fitting for each channel over whole partition from cal data
- combining of all channels into single pars files with associated plot and results files
- running build hit over all channels using par file
"""

from scripts.util.pars_loading import pars_catalog
from scripts.util.create_pars_keylist import pars_key_resolve
from scripts.util.utils import filelist_path, par_pht_path, set_last_rule_name
from scripts.util.patterns import (
    get_pattern_pars_tmp_channel,
    get_pattern_plts_tmp_channel,
    get_pattern_log_channel,
    get_pattern_par_pht,
    get_pattern_plts,
    get_pattern_tier,
    get_pattern_pars_tmp,
    get_pattern_log,
    get_pattern_pars,
)

ds.pars_key_resolve.write_par_catalog(
    ["-*-*-*-cal"],
    os.path.join(pars_path(setup), "pht", "validity.jsonl"),
    get_pattern_tier_raw(setup),
    {"cal": ["par_pht"], "lar": ["par_pht"]},
)

intier = "dsp"


rule pht_checkpoint:
    input:
        files=lambda wildcards: read_filelist_cal(wildcards, intier),
    output:
        get_pattern_pars_tmp_channel(setup, "pht", "check"),
    shell:
        "touch {output}"


qc_pht_rules = {}
for key, dataset in part.datasets.items():
    for partition in dataset.keys():

        rule:
            input:
                cal_files=part.get_filelists(partition, key, intier),
                fft_files=part.get_filelists(partition, key, intier, datatype="fft"),
                pulser_files=[
                    file.replace("pht", "tcm")
                    for file in part.get_par_files(
                        f"{par_pht_path(setup)}/validity.jsonl",
                        partition,
                        key,
                        tier="pht",
                        name="pulser_ids",
                    )
                ],
                check_files=part.get_par_files(
                    f"{par_pht_path(setup)}/validity.jsonl",
                    partition,
                    key,
                    tier="pht",
                    name="check",
                ),
            wildcard_constraints:
                channel=part.get_wildcard_constraints(partition, key),
            params:
                datatype="cal",
                channel="{channel}" if key == "default" else key,
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
                        name="qc",
                    )
                ],
                plot_file=[
                    temp(file)
                    for file in part.get_plt_files(
                        f"{par_pht_path(setup)}/validity.jsonl",
                        partition,
                        key,
                        tier="pht",
                        name="qc",
                    )
                ],
            log:
                part.get_log_file(
                    f"{par_pht_path(setup)}/validity.jsonl",
                    partition,
                    key,
                    "pht",
                    name="par_pht_qc",
                ),
            group:
                "par-pht"
            resources:
                mem_swap=len(part.get_filelists(partition, key, intier)) * 20,
                runtime=300,
            shell:
                "{swenv} python3 -B "
                f"{basedir}/../scripts/pars_pht_qc.py "
                "--log {log} "
                "--configs {configs} "
                "--datatype {params.datatype} "
                "--timestamp {params.timestamp} "
                "--channel {params.channel} "
                "--save_path {output.hit_pars} "
                "--plot_path {output.plot_file} "
                "--pulser_files {input.pulser_files} "
                "--fft_files {input.fft_files} "
                "--cal_files {input.cal_files}"

        set_last_rule_name(workflow, f"{key}-{partition}-build_pht_qc")

        if key in qc_pht_rules:
            qc_pht_rules[key].append(list(workflow.rules)[-1])
        else:
            qc_pht_rules[key] = [list(workflow.rules)[-1]]


# Merged energy and a/e supercalibrations to reduce number of rules as they have same inputs/outputs
# This rule builds the a/e calibration using the calibration dsp files for the whole partition
rule build_pht_qc:
    input:
        cal_files=os.path.join(
            filelist_path(setup),
            "all-{experiment}-{period}-{run}-cal-" + f"{intier}.filelist",
        ),
        fft_files=os.path.join(
            filelist_path(setup),
            "all-{experiment}-{period}-{run}-fft-" + f"{intier}.filelist",
        ),
        pulser_file=get_pattern_pars_tmp_channel(setup, "tcm", "pulser_ids"),
        check_file=get_pattern_pars_tmp_channel(setup, "pht", "check"),
    params:
        datatype="cal",
        channel="{channel}",
        timestamp="{timestamp}",
    output:
        hit_pars=temp(get_pattern_pars_tmp_channel(setup, "pht", "qc")),
        plot_file=temp(get_pattern_plts_tmp_channel(setup, "pht", "qc")),
    log:
        get_pattern_log_channel(setup, "pars_pht_qc"),
    group:
        "par-pht"
    resources:
        mem_swap=60,
        runtime=300,
    shell:
        "{swenv} python3 -B "
        f"{basedir}/../scripts/pars_pht_qc.py "
        "--log {log} "
        "--configs {configs} "
        "--datatype {params.datatype} "
        "--timestamp {params.timestamp} "
        "--channel {params.channel} "
        "--save_path {output.hit_pars} "
        "--plot_path {output.plot_file} "
        "--pulser_files {input.pulser_files} "
        "--fft_files {input.fft_files} "
        "--cal_files {input.cal_files}"


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
        files=os.path.join(
            filelist_path(setup),
            "all-{experiment}-{period}-{run}-cal-" + f"{intier}.filelist",
        ),
        pulser=get_pattern_pars_tmp_channel(setup, "tcm", "pulser_ids"),
        pht_dict=get_pattern_pars_tmp_channel(setup, "pht", "qc"),
        inplots=get_pattern_plts_tmp_channel(setup, "pht", "qc"),
        ctc_dict=ancient(
            lambda wildcards: pars_catalog.get_par_file(
                setup, wildcards.timestamp, intier
            )
        ),
    params:
        timestamp="{timestamp}",
        datatype="cal",
        channel="{channel}",
        tier="pht",
    output:
        ecal_file=temp(get_pattern_pars_tmp_channel(setup, "pht", "energy_cal")),
        results_file=temp(
            get_pattern_pars_tmp_channel(
                setup, "pht", "energy_cal_objects", extension="pkl"
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
        "{swenv} python3 -B "
        f"{workflow.source_path('../scripts/pars_hit_ecal.py')} "
        "--log {log} "
        "--datatype {params.datatype} "
        "--timestamp {params.timestamp} "
        "--channel {params.channel} "
        "--configs {configs} "
        "--tier {params.tier} "
        "--metadata {meta} "
        "--plot_path {output.plot_file} "
        "--results_path {output.results_file} "
        "--save_path {output.ecal_file} "
        "--inplot_dict {input.inplots} "
        "--in_hit_dict {input.pht_dict} "
        "--ctc_dict {input.ctc_dict} "
        "--pulser_file {input.pulser} "
        "--files {input.files}"


part_pht_rules = {}
for key, dataset in part.datasets.items():
    for partition in dataset.keys():

        rule:
            input:
                files=part.get_filelists(partition, key, intier),
                pulser_files=[
                    file.replace("pht", "tcm")
                    for file in part.get_par_files(
                        f"{par_pht_path(setup)}/validity.jsonl",
                        partition,
                        key,
                        tier="pht",
                        name="pulser_ids",
                    )
                ],
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
                    name="energy_cal_objects",
                    extension="pkl",
                ),
                inplots=part.get_plt_files(
                    f"{par_pht_path(setup)}/validity.jsonl",
                    partition,
                    key,
                    tier="pht",
                    name="energy_cal",
                ),
            wildcard_constraints:
                channel=part.get_wildcard_constraints(partition, key),
            params:
                datatype="cal",
                channel="{channel}" if key == "default" else key,
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
                        name="partcal",
                    )
                ],
                partcal_results=[
                    temp(file)
                    for file in part.get_par_files(
                        f"{par_pht_path(setup)}/validity.jsonl",
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
                        f"{par_pht_path(setup)}/validity.jsonl",
                        partition,
                        key,
                        tier="pht",
                        name="partcal",
                    )
                ],
            log:
                part.get_log_file(
                    f"{par_pht_path(setup)}/validity.jsonl",
                    partition,
                    key,
                    "pht",
                    name="par_pht_partcal",
                ),
            group:
                "par-pht"
            resources:
                mem_swap=len(part.get_filelists(partition, key, intier)) * 20,
                runtime=300,
            shell:
                "{swenv} python3 -B "
                f"{basedir}/../scripts/pars_pht_partcal.py "
                "--log {log} "
                "--configs {configs} "
                "--datatype {params.datatype} "
                "--timestamp {params.timestamp} "
                "--inplots {input.inplots} "
                "--channel {params.channel} "
                "--metadata {meta} "
                "--fit_results {output.partcal_results} "
                "--eres_file {input.eres_file} "
                "--hit_pars {output.hit_pars} "
                "--plot_file {output.plot_file} "
                "--ecal_file {input.ecal_file} "
                "--pulser_files {input.pulser_files} "
                "--input_files {input.files}"

        set_last_rule_name(
            workflow, f"{key}-{partition}-build_pht_energy_super_calibrations"
        )

        if key in part_pht_rules:
            part_pht_rules[key].append(list(workflow.rules)[-1])
        else:
            part_pht_rules[key] = [list(workflow.rules)[-1]]


# Merged energy and a/e supercalibrations to reduce number of rules as they have same inputs/outputs
# This rule builds the a/e calibration using the calibration dsp files for the whole partition
rule build_pht_energy_super_calibrations:
    input:
        files=os.path.join(
            filelist_path(setup),
            "all-{experiment}-{period}-{run}-cal" + f"-{intier}.filelist",
        ),
        pulser_files=get_pattern_pars_tmp_channel(setup, "tcm", "pulser_ids"),
        ecal_file=get_pattern_pars_tmp_channel(setup, "pht", "energy_cal"),
        eres_file=get_pattern_pars_tmp_channel(
            setup, "pht", "energy_cal_objects", extension="pkl"
        ),
        inplots=get_pattern_plts_tmp_channel(setup, "pht", "energy_cal"),
    params:
        datatype="cal",
        channel="{channel}",
        timestamp="{timestamp}",
    output:
        hit_pars=temp(get_pattern_pars_tmp_channel(setup, "pht", "partcal")),
        partcal_results=temp(
            get_pattern_pars_tmp_channel(
                setup, "pht", "partcal_objects", extension="pkl"
            )
        ),
        plot_file=temp(get_pattern_plts_tmp_channel(setup, "pht", "partcal")),
    log:
        get_pattern_log_channel(setup, "pars_pht_partcal"),
    group:
        "par-pht"
    resources:
        mem_swap=60,
        runtime=300,
    shell:
        "{swenv} python3 -B "
        f"{basedir}/../scripts/pars_pht_partcal.py "
        "--log {log} "
        "--configs {configs} "
        "--datatype {params.datatype} "
        "--timestamp {params.timestamp} "
        "--channel {params.channel} "
        "--metadata {meta} "
        "--inplots {input.inplots} "
        "--fit_results {output.partcal_results} "
        "--eres_file {input.eres_file} "
        "--hit_pars {output.hit_pars} "
        "--plot_file {output.plot_file} "
        "--ecal_file {input.ecal_file} "
        "--pulser_files {input.pulser_files} "
        "--input_files {input.files}"


fallback_pht_rule = list(workflow.rules)[-1]

rule_order_list = []
ordered = OrderedDict(part_pht_rules)
ordered.move_to_end("default")
for key, items in ordered.items():
    rule_order_list += [item.name for item in items]
rule_order_list.append(fallback_pht_rule.name)
workflow._ruleorder.add(*rule_order_list)  # [::-1]

part_pht_rules = {}
for key, dataset in part.datasets.items():
    for partition in dataset.keys():

        rule:
            input:
                files=part.get_filelists(partition, key, intier),
                pulser_files=[
                    file.replace("pht", "tcm")
                    for file in part.get_par_files(
                        f"{par_pht_path(setup)}/validity.jsonl",
                        partition,
                        key,
                        tier="pht",
                        name="pulser_ids",
                    )
                ],
                ecal_file=part.get_par_files(
                    f"{par_pht_path(setup)}/validity.jsonl",
                    partition,
                    key,
                    tier="pht",
                    name="partcal",
                ),
                eres_file=part.get_par_files(
                    f"{par_pht_path(setup)}/validity.jsonl",
                    partition,
                    key,
                    tier="pht",
                    name="partcal_objects",
                    extension="pkl",
                ),
                inplots=part.get_plt_files(
                    f"{par_pht_path(setup)}/validity.jsonl",
                    partition,
                    key,
                    tier="pht",
                    name="partcal",
                ),
            wildcard_constraints:
                channel=part.get_wildcard_constraints(partition, key),
            params:
                datatype="cal",
                channel="{channel}" if key == "default" else key,
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
                        name="aoecal",
                    )
                ],
                aoe_results=[
                    temp(file)
                    for file in part.get_par_files(
                        f"{par_pht_path(setup)}/validity.jsonl",
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
                        f"{par_pht_path(setup)}/validity.jsonl",
                        partition,
                        key,
                        tier="pht",
                        name="aoecal",
                    )
                ],
            log:
                part.get_log_file(
                    f"{par_pht_path(setup)}/validity.jsonl",
                    partition,
                    key,
                    "pht",
                    name="par_pht_aoe",
                ),
            group:
                "par-pht"
            resources:
                mem_swap=len(part.get_filelists(partition, key, intier)) * 20,
                runtime=300,
            shell:
                "{swenv} python3 -B "
                f"{basedir}/../scripts/pars_pht_aoecal.py "
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
                "--pulser_files {input.pulser_files} "
                "--input_files {input.files}"

        set_last_rule_name(
            workflow, f"{key}-{partition}-build_pht_aoe_calibrations"
        )

        if key in part_pht_rules:
            part_pht_rules[key].append(list(workflow.rules)[-1])
        else:
            part_pht_rules[key] = [list(workflow.rules)[-1]]


# Merged energy and a/e supercalibrations to reduce number of rules as they have same inputs/outputs
# This rule builds the a/e calibration using the calibration dsp files for the whole partition
rule build_pht_aoe_calibrations:
    input:
        files=os.path.join(
            filelist_path(setup),
            "all-{experiment}-{period}-{run}-cal-" + f"{intier}.filelist",
        ),
        pulser_files=get_pattern_pars_tmp_channel(setup, "tcm", "pulser_ids"),
        ecal_file=get_pattern_pars_tmp_channel(setup, "pht", "partcal"),
        eres_file=get_pattern_pars_tmp_channel(
            setup, "pht", "partcal_objects", extension="pkl"
        ),
        inplots=get_pattern_plts_tmp_channel(setup, "pht", "partcal"),
    params:
        datatype="cal",
        channel="{channel}",
        timestamp="{timestamp}",
    output:
        hit_pars=temp(get_pattern_pars_tmp_channel(setup, "pht", "aoecal")),
        aoe_results=temp(
            get_pattern_pars_tmp_channel(
                setup, "pht", "aoecal_objects", extension="pkl"
            )
        ),
        plot_file=temp(get_pattern_plts_tmp_channel(setup, "pht", "aoecal")),
    log:
        get_pattern_log_channel(setup, "pars_pht_aoe_cal"),
    group:
        "par-pht"
    resources:
        mem_swap=60,
        runtime=300,
    shell:
        "{swenv} python3 -B "
        f"{basedir}/../scripts/pars_pht_aoecal.py "
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
        "--pulser_files {input.pulser_files} "
        "--input_files {input.files}"


fallback_pht_rule = list(workflow.rules)[-1]

rule_order_list = []
ordered = OrderedDict(part_pht_rules)
ordered.move_to_end("default")
for key, items in ordered.items():
    rule_order_list += [item.name for item in items]
rule_order_list.append(fallback_pht_rule.name)
workflow._ruleorder.add(*rule_order_list)  # [::-1]

part_pht_rules = {}
for key, dataset in part.datasets.items():
    for partition in dataset.keys():

        rule:
            input:
                files=part.get_filelists(partition, key, intier),
                pulser_files=[
                    file.replace("pht", "tcm")
                    for file in part.get_par_files(
                        f"{par_pht_path(setup)}/validity.jsonl",
                        partition,
                        key,
                        tier="pht",
                        name="pulser_ids",
                    )
                ],
                ecal_file=part.get_par_files(
                    f"{par_pht_path(setup)}/validity.jsonl",
                    partition,
                    key,
                    tier="pht",
                    name="aoecal",
                ),
                eres_file=part.get_par_files(
                    f"{par_pht_path(setup)}/validity.jsonl",
                    partition,
                    key,
                    tier="pht",
                    name="aoecal_objects",
                    extension="pkl",
                ),
                inplots=part.get_plt_files(
                    f"{par_pht_path(setup)}/validity.jsonl",
                    partition,
                    key,
                    tier="pht",
                    name="aoecal",
                ),
            wildcard_constraints:
                channel=part.get_wildcard_constraints(partition, key),
            params:
                datatype="cal",
                channel="{channel}" if key == "default" else key,
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
                lq_results=[
                    temp(file)
                    for file in part.get_par_files(
                        f"{par_pht_path(setup)}/validity.jsonl",
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
                    name="par_pht_lq",
                ),
            group:
                "par-pht"
            resources:
                mem_swap=len(part.get_filelists(partition, key, intier)) * 20,
                runtime=300,
            shell:
                "{swenv} python3 -B "
                f"{basedir}/../scripts/pars_pht_lqcal.py "
                "--log {log} "
                "--configs {configs} "
                "--datatype {params.datatype} "
                "--timestamp {params.timestamp} "
                "--inplots {input.inplots} "
                "--channel {params.channel} "
                "--lq_results {output.lq_results} "
                "--eres_file {input.eres_file} "
                "--hit_pars {output.hit_pars} "
                "--plot_file {output.plot_file} "
                "--ecal_file {input.ecal_file} "
                "--pulser_files {input.pulser_files} "
                "--input_files {input.files}"

        set_last_rule_name(workflow, f"{key}-{partition}-build_pht_lq_calibration")

        if key in part_pht_rules:
            part_pht_rules[key].append(list(workflow.rules)[-1])
        else:
            part_pht_rules[key] = [list(workflow.rules)[-1]]


# This rule builds the lq calibration using the calibration dsp files for the whole partition
rule build_pht_lq_calibration:
    input:
        files=os.path.join(
            filelist_path(setup),
            "all-{experiment}-{period}-{run}-cal-" + f"{intier}.filelist",
        ),
        pulser_files=get_pattern_pars_tmp_channel(setup, "tcm", "pulser_ids"),
        ecal_file=get_pattern_pars_tmp_channel(setup, "pht", "aoecal"),
        eres_file=get_pattern_pars_tmp_channel(
            setup, "pht", "aoecal_objects", extension="pkl"
        ),
        inplots=get_pattern_plts_tmp_channel(setup, "pht", "aoecal"),
    params:
        datatype="cal",
        channel="{channel}",
        timestamp="{timestamp}",
    output:
        hit_pars=temp(get_pattern_pars_tmp_channel(setup, "pht")),
        lq_results=temp(
            get_pattern_pars_tmp_channel(setup, "pht", "objects", extension="pkl")
        ),
        plot_file=temp(get_pattern_plts_tmp_channel(setup, "pht")),
    log:
        get_pattern_log_channel(setup, "pars_pht_lq_cal"),
    group:
        "par-pht"
    resources:
        mem_swap=60,
        runtime=300,
    shell:
        "{swenv} python3 -B "
        f"{basedir}/../scripts/pars_pht_lqcal.py "
        "--log {log} "
        "--configs {configs} "
        "--datatype {params.datatype} "
        "--timestamp {params.timestamp} "
        "--inplots {input.inplots} "
        "--channel {params.channel} "
        "--lq_results {output.lq_results} "
        "--eres_file {input.eres_file} "
        "--hit_pars {output.hit_pars} "
        "--plot_file {output.plot_file} "
        "--ecal_file {input.ecal_file} "
        "--pulser_files {input.pulser_files} "
        "--input_files {input.files}"


fallback_pht_rule = list(workflow.rules)[-1]

rule_order_list = []
ordered = OrderedDict(part_pht_rules)
ordered.move_to_end("default")
for key, items in ordered.items():
    rule_order_list += [item.name for item in items]
rule_order_list.append(fallback_pht_rule.name)
workflow._ruleorder.add(*rule_order_list)  # [::-1]


rule build_pars_pht_objects:
    input:
        lambda wildcards: read_filelist_pars_cal_channel(
            wildcards,
            "pht_objects_pkl",
        ),
    output:
        get_pattern_pars(
            setup,
            "pht",
            name="objects",
            extension="dir",
            check_in_cycle=check_in_cycle,
        ),
    group:
        "merge-hit"
    shell:
        "{swenv} python3 -B "
        f"{basedir}/../scripts/merge_channels.py "
        "--input {input} "
        "--output {output} "


rule build_plts_pht:
    input:
        lambda wildcards: read_filelist_plts_cal_channel(wildcards, "pht"),
    output:
        get_pattern_plts(setup, "pht"),
    group:
        "merge-hit"
    shell:
        "{swenv} python3 -B "
        f"{basedir}/../scripts/merge_channels.py "
        "--input {input} "
        "--output {output} "


rule build_pars_pht:
    input:
        infiles=lambda wildcards: read_filelist_pars_cal_channel(wildcards, "pht"),
        plts=get_pattern_plts(setup, "pht"),
        objects=get_pattern_pars(
            setup,
            "pht",
            name="objects",
            extension="dir",
            check_in_cycle=check_in_cycle,
        ),
    output:
        get_pattern_pars(setup, "pht", check_in_cycle=check_in_cycle),
    group:
        "merge-hit"
    shell:
        "{swenv} python3 -B "
        f"{basedir}/../scripts/merge_channels.py "
        "--input {input.infiles} "
        "--output {output} "


rule build_pht:
    input:
        dsp_file=get_pattern_tier(setup, intier, check_in_cycle=False),
        pars_file=lambda wildcards: pars_catalog.get_par_file(
            setup, wildcards.timestamp, "pht"
        ),
    output:
        tier_file=get_pattern_tier(setup, "pht", check_in_cycle=check_in_cycle),
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
