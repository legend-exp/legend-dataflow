from scripts.util.pars_loading import ParsCatalog
from scripts.util.create_pars_keylist import ParsKeyResolve
from scripts.util.utils import filelist_path, set_last_rule_name
from scripts.util.patterns import (
    get_pattern_pars_tmp_channel,
    get_pattern_plts_tmp_channel,
    get_pattern_log_channel,
    get_pattern_plts,
    get_pattern_tier,
    get_pattern_pars_tmp,
    get_pattern_log,
    get_pattern_pars,
)


pht_fast_rules = {}
for key, dataset in part.datasets.items():
    for partition in dataset.keys():

        rule:
            input:
                files=part.get_filelists(partition, key, intier),
                pulser_files=[
                    file.replace("par_pht", "par_tcm")
                    for file in part.get_par_files(
                        pht_par_catalog,
                        partition,
                        key,
                        tier="pht",
                        name="pulser_ids",
                    )
                ],
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
                channel=part.get_wildcard_constraints(partition, key),
            params:
                datatype="cal",
                channel="{channel}" if key == "default" else key,
                timestamp=part.get_timestamp(
                    pht_par_catalog, partition, key, tier="pht"
                ),
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
                partcal_results=[
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
                    name="par_pht_fast",
                ),
            group:
                "par-pht"
            resources:
                mem_swap=len(part.get_filelists(partition, key, intier)) * 12,
                runtime=300,
            shell:
                "{swenv} python3 -B "
                f"{basedir}/../scripts/pars_pht_fast.py "
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

        set_last_rule_name(workflow, f"{key}-{partition}-par_pht_fast")
        slow_rule = workflow._rules[f"{key}-{partition}-build_pht_lq_calibration"]

        if key in pht_fast_rules:
            pht_fast_rules[key] += [list(workflow.rules)[-1], slow_rule]
        else:
            pht_fast_rules[key] = [list(workflow.rules)[-1], slow_rule]


# Merged energy and a/e supercalibrations to reduce number of rules as they have same inputs/outputs
# This rule builds the a/e calibration using the calibration dsp files for the whole partition
rule par_pht_fast:
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
        hit_pars=temp(get_pattern_pars_tmp_channel(setup, "pht")),
        partcal_results=temp(
            get_pattern_pars_tmp_channel(setup, "pht", "objects", extension="pkl")
        ),
        plot_file=temp(get_pattern_plts_tmp_channel(setup, "pht")),
    log:
        get_pattern_log_channel(setup, "par_pht_fast"),
    group:
        "par-pht"
    resources:
        mem_swap=50,
        runtime=300,
    shell:
        "{swenv} python3 -B "
        "{basedir}/../scripts/pars_pht_fast.py "
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
slow_fallback = workflow._rules["build_pht_lq_calibration"]

rule_order_list = []
ordered = OrderedDict(pht_fast_rules)
ordered.move_to_end("default")
for key, items in ordered.items():
    rule_order_list += [item.name for item in items]
rule_order_list.append(fallback_pht_rule.name)
rule_order_list.append(slow_fallback.name)
workflow._ruleorder.add(*rule_order_list)  # [::-1]
