"""
Snakemake rules for processing evt tier.
"""

from legenddataflow.pars_loading import ParsCatalog
from legenddataflow.patterns import (
    get_pattern_tier,
    get_pattern_log,
    get_pattern_pars,
    get_pattern_log_concat,
)
from legenddataflow.execenv import execenv_pyexe


rule build_evt:
    input:
        dsp_file=get_pattern_tier(config, "dsp", check_in_cycle=False),
        hit_file=get_pattern_tier(config, "hit", check_in_cycle=False),
        tcm_file=get_pattern_tier(config, "tcm", check_in_cycle=False),
        ann_file=lambda wildcards: (
            []
            if int(wildcards["period"][1:]) > 11
            else get_pattern_tier(config, "ann", check_in_cycle=False)
        ),
        par_files=lambda wildcards: hit_par_catalog.get_par_file(
            config, wildcards.timestamp, "hit"
        ),
        xtalk_matrix=lambda wildcards: get_input_par_file(
            tier="evt", wildcards=wildcards, name="xtc"
        ),
    output:
        get_pattern_tier(config, "evt", check_in_cycle=check_in_cycle),
    params:
        timestamp="{timestamp}",
        datatype="{datatype}",
        tier="evt",
        ro_input=lambda _, input: {k: ro(v) for k, v in input.items()},
    log:
        get_pattern_log(config, f"tier_evt", time),
    group:
        "tier-evt"
    resources:
        runtime=300,
        mem_swap=50,
    run:
        shell_string = (
            execenv_pyexe(config, "build-tier-evt") + f"--configs {ro(configs)} "
            f"--metadata {ro(meta)} "
            "--log {log} "
            "--tier {params.tier} "
            "--datatype {params.datatype} "
            "--timestamp {params.timestamp} "
            "--xtc-file {params.ro_input[xtalk_matrix]} "
            "--par-files {params.ro_input[par_files]} "
            "--hit-file {params.ro_input[hit_file]} "
            "--tcm-file {params.ro_input[tcm_file]} "
            "--dsp-file {params.ro_input[dsp_file]} "
            "--output {output} "
        )
        if input.ann_file is not None:
            shell_string += "--ann-file {params.ro_input[ann_file]} "

        shell(shell_string)


rule build_pet:
    input:
        dsp_file=get_pattern_tier(config, "psp", check_in_cycle=False),
        hit_file=get_pattern_tier(config, "pht", check_in_cycle=False),
        tcm_file=get_pattern_tier(config, "tcm", check_in_cycle=False),
        ann_file=lambda wildcards: (
            []
            if int(wildcards["period"][1:]) > 11
            else get_pattern_tier(config, "pan", check_in_cycle=False)
        ),
        par_files=lambda wildcards: pht_par_catalog.get_par_file(
            config, wildcards.timestamp, "pht"
        ),
        xtalk_matrix=lambda wildcards: get_input_par_file(
            tier="pet", wildcards=wildcards, name="xtc"
        ),
    output:
        get_pattern_tier(config, "pet", check_in_cycle=check_in_cycle),
    params:
        timestamp="{timestamp}",
        datatype="{datatype}",
        tier="pet",
        ro_input=lambda _, input: {k: ro(v) for k, v in input.items()},
    log:
        get_pattern_log(config, f"tier_pet", time),
    group:
        "tier-evt"
    resources:
        runtime=300,
        mem_swap=50,
    run:
        shell_string = (
            execenv_pyexe(config, "build-tier-evt") + f"--configs {ro(configs)} "
            f"--metadata {ro(meta)} "
            "--log {log} "
            "--tier {params.tier} "
            "--datatype {params.datatype} "
            "--timestamp {params.timestamp} "
            "--xtc-file {params.ro_input[xtalk_matrix]} "
            "--par-files {params.ro_input[par_files]} "
            "--hit-file {params.ro_input[hit_file]} "
            "--tcm-file {params.ro_input[tcm_file]} "
            "--dsp-file {params.ro_input[dsp_file]} "
            "--output {output} "
        )
        if input.ann_file is not None:
            shell_string += "--ann-file {params.ro_input[ann_file]} "

        shell(shell_string)


for evt_tier in ("evt", "pet"):

    rule:
        wildcard_constraints:
            timestamp=r"(?!\d{8}T\d{6}Z)",
        input:
            lambda wildcards: sorted(
                get_filelist_full_wildcards(
                    wildcards,
                    config,
                    get_pattern_tier_raw(config),
                    tier,
                    ignore_keys_file=os.path.join(configs, "ignore_keys.keylist"),
                )
            ),
        output:
            get_pattern_tier(
                config, f"{evt_tier}_concat", check_in_cycle=check_in_cycle
            ),
        params:
            timestamp="all",
            datatype="{datatype}",
            ro_input=lambda _, input: utils.as_ro(config, input),
        log:
            get_pattern_log_concat(config, f"tier_{evt_tier}_concat", time),
        group:
            "tier-evt"
        shell:
            execenv_pyexe(config, "lh5concat") + "--verbose --overwrite "
            "--output {output} "
            "-- {params.ro_input} &> {log}"

    set_last_rule_name(workflow, f"concat_{evt_tier}")
