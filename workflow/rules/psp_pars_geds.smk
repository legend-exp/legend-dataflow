"""
Snakemake rules for processing psp (partition dsp) tier data.
- extraction of calibration curves(s) for each run for each channel from cal data
- extraction of psd calibration parameters and partition level energy fitting for each channel over whole partition from cal data
"""

from legenddataflow.utils import set_last_rule_name
from legenddataflow.create_pars_keylist import ParsKeyResolve
from legenddataflow.patterns import (
    get_pattern_pars_tmp_channel,
    get_pattern_plts_tmp_channel,
    get_pattern_log_channel,
    get_pattern_log,
    get_pattern_pars,
    get_pattern_tier,
)
from legenddataflow.execenv import execenv_smk_py_script

psp_par_catalog = ParsKeyResolve.get_par_catalog(
    ["-*-*-*-cal"],
    get_pattern_tier(config, "raw", check_in_cycle=False),
    {"cal": ["par_psp"], "lar": ["par_psp"]},
)


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
            wildcard_constraints:
                channel=part.get_wildcard_constraints(partition, key),
            params:
                datatype="cal",
                channel="{channel}" if key == "default" else key,
                timestamp=part.get_timestamp(
                    psp_par_catalog, partition, key, tier="psp"
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
            group:
                "par-psp"
            resources:
                runtime=300,
            shell:
                f'{execenv_smk_py_script(config, "par-geds-psp-average")}'
                "--log {log} "
                "--configs {configs} "
                "--datatype {params.datatype} "
                "--timestamp {params.timestamp} "
                "--channel {params.channel} "
                "--in_plots {input.dsp_plots} "
                "--out_plots {output.psp_plots} "
                "--in_obj {input.dsp_objs} "
                "--out_obj {output.psp_objs} "
                "--input {input.dsp_pars} "
                "--output {output.psp_pars} "

        set_last_rule_name(workflow, f"{key}-{partition}-build_par_psp")

        if key in psp_rules:
            psp_rules[key].append(list(workflow.rules)[-1])
        else:
            psp_rules[key] = [list(workflow.rules)[-1]]


# Merged energy and a/e supercalibrations to reduce number of rules as they have same inputs/outputs
# This rule builds the a/e calibration using the calibration dsp files for the whole partition
rule build_par_psp_fallback:
    input:
        dsp_pars=get_pattern_pars_tmp_channel(config, "dsp", "eopt"),
        dsp_objs=get_pattern_pars_tmp_channel(config, "dsp", "objects", extension="pkl"),
        dsp_plots=get_pattern_plts_tmp_channel(config, "dsp"),
    params:
        datatype="cal",
        channel="{channel}",
        timestamp="{timestamp}",
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
    shell:
        f'{execenv_smk_py_script(config, "par-geds-psp-average")}'
        "--log {log} "
        "--configs {configs} "
        "--datatype {params.datatype} "
        "--timestamp {params.timestamp} "
        "--channel {params.channel} "
        "--in_plots {input.dsp_plots} "
        "--out_plots {output.psp_plots} "
        "--in_obj {input.dsp_objs} "
        "--out_obj {output.psp_objs} "
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


rule build_svm_psp:
    input:
        hyperpars=lambda wildcards: get_input_par_file(
            wildcards, "psp", "svm_hyperpars"
        ),
        train_data=lambda wildcards: str(
            get_input_par_file(wildcards, "psp", "svm_hyperpars")
        ).replace("hyperpars.json", "train.lh5"),
    output:
        dsp_pars=get_pattern_pars(config, "psp", "svm", "pkl"),
    log:
        str(get_pattern_log(config, "pars_psp_svm", time)).replace("{datatype}", "cal"),
    group:
        "par-dsp-svm"
    resources:
        runtime=300,
    shell:
        f'{execenv_smk_py_script(config, "par-geds-dsp-svm-build")}'
        "--log {log} "
        "--train_data {input.train_data} "
        "--train_hyperpars {input.hyperpars} "
        "--output_file {output.dsp_pars}"


rule build_pars_psp_svm:
    input:
        dsp_pars=get_pattern_pars_tmp_channel(config, "psp_eopt"),
        svm_model=get_pattern_pars(config, "psp", "svm", "pkl"),
    output:
        dsp_pars=temp(get_pattern_pars_tmp_channel(config, "psp")),
    log:
        get_pattern_log_channel(config, "pars_dsp_svm", time),
    group:
        "par-dsp"
    resources:
        runtime=300,
    shell:
        f'{execenv_smk_py_script(config, "par-geds-dsp-svm")}'
        "--log {log} "
        "--input_file {input.dsp_pars} "
        "--output_file {output.dsp_pars} "
        "--svm_file {input.svm_model}"
