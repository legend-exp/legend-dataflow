from legenddataflow.pars_loading import ParsCatalog
from legenddataflow.create_pars_keylist import ParsKeyResolve
from legenddataflow.utils import set_last_rule_name
from legenddataflow.paths import filelist_path
from legenddataflow.patterns import (
    get_pattern_pars_tmp_channel,
    get_pattern_plts_tmp_channel,
    get_pattern_log_channel,
    get_pattern_plts,
    get_pattern_tier,
    get_pattern_pars_tmp,
    get_pattern_log,
    get_pattern_pars,
)
from legenddataflow.execenv import execenv_pyexe

intier = "psp"


qc_pht_rules = {}
for key, dataset in part.datasets.items():
    for partition in dataset.keys():

        rule:
            input:
                phy_files=part.get_filelists(partition, key, intier, datatype="phy"),
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
                        name="qcphy",
                    )
                ],
                plot_file=[
                    temp(file)
                    for file in part.get_plt_files(
                        pht_par_catalog,
                        partition,
                        key,
                        tier="pht",
                        name="qcphy",
                    )
                ],
            log:
                part.get_log_file(
                    pht_par_catalog,
                    partition,
                    key,
                    "pht",
                    time,
                    name="par_pht_qc_phy",
                ),
            group:
                "par-pht"
            resources:
                mem_swap=len(part.get_filelists(partition, key, intier)) * 20,
                runtime=300,
            shell:
                execenv_pyexe(config, "par-geds-pht-qc-phy") + "--log {log} "
                "--configs {configs} "
                "--datatype {params.datatype} "
                "--timestamp {params.timestamp} "
                "--channel {params.channel} "
                "--save-path {output.hit_pars} "
                "--plot-path {output.plot_file} "
                "--phy-files {input.phy_files}"

        set_last_rule_name(workflow, f"{key}-{partition}-build_pht_qc_phy")

        if key in qc_pht_rules:
            qc_pht_rules[key].append(list(workflow.rules)[-1])
        else:
            qc_pht_rules[key] = [list(workflow.rules)[-1]]


# Merged energy and a/e supercalibrations to reduce number of rules as they have same inputs/outputs
# This rule builds the a/e calibration using the calibration dsp files for the whole partition
rule build_pht_qc_phy:
    input:
        phy_files=os.path.join(
            filelist_path(config),
            "all-{experiment}-{period}-{run}-phy-" + f"{intier}.filelist",
        ),
    params:
        datatype="cal",
        channel="{channel}",
        timestamp="{timestamp}",
    output:
        hit_pars=temp(get_pattern_pars_tmp_channel(config, "pht", "qcphy")),
        plot_file=temp(get_pattern_plts_tmp_channel(config, "pht", "qcphy")),
    log:
        get_pattern_log_channel(config, "pars_pht_qc_phy", time),
    group:
        "par-pht"
    resources:
        mem_swap=60,
        runtime=300,
    shell:
        execenv_pyexe(config, "par-geds-pht-qc-phy") + "--log {log} "
        "--configs {configs} "
        "--datatype {params.datatype} "
        "--timestamp {params.timestamp} "
        "--channel {params.channel} "
        "--save-path {output.hit_pars} "
        "--plot-path {output.plot_file} "
        "--phy-files {input.phy_files}"


fallback_qc_rule = list(workflow.rules)[-1]

rule_order_list = []
ordered = OrderedDict(qc_pht_rules)
ordered.move_to_end("default")
for key, items in ordered.items():
    rule_order_list += [item.name for item in items]
rule_order_list.append(fallback_qc_rule.name)
workflow._ruleorder.add(*rule_order_list)  # [::-1]


rule build_plts_pht_phy:
    input:
        lambda wildcards: get_plt_chanlist(
            config,
            f"all-{wildcards.experiment}-{wildcards.period}-{wildcards.run}-cal-{wildcards.timestamp}-channels",
            "pht",
            basedir,
            det_status,
            chan_maps,
            name="qcphy",
        ),
    output:
        get_pattern_plts(config, "pht", "qc_phy"),
    group:
        "merge-hit"
    shell:
        execenv_pyexe(config, "merge-channels") + "--input {input} "
        "--output {output} "


rule build_pars_pht_phy:
    input:
        infiles=lambda wildcards: get_par_chanlist(
            config,
            f"all-{wildcards.experiment}-{wildcards.period}-{wildcards.run}-cal-{wildcards.timestamp}-channels",
            "pht",
            basedir,
            det_status,
            chan_maps,
            name="qcphy",
        ),
        plts=get_pattern_plts(config, "pht", "qc_phy"),
    output:
        get_pattern_pars(config, "pht", name="qc_phy", check_in_cycle=check_in_cycle),
    group:
        "merge-hit"
    shell:
        execenv_pyexe(config, "merge-channels") + "--input {input.infiles} "
        "--output {output} "
