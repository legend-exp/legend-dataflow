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


rule build_evt:
    input:
        dsp_file=get_pattern_tier(setup, "dsp", check_in_cycle=False),
        hit_file=get_pattern_tier(setup, "hit", check_in_cycle=False),
        tcm_file=get_pattern_tier(setup, "tcm", check_in_cycle=False),
        ann_file=lambda wildcards: (
            None
            if int(wildcards["period"][1:]) > 11
            else get_pattern_tier(setup, "ann", check_in_cycle=False)
        ),
        par_files=lambda wildcards: ParsCatalog.get_par_file(
            setup, wildcards.timestamp, "hit"
        ),
        xtalk_matrix=lambda wildcards: get_input_par_file(
            tier="evt", wildcards=wildcards, name="xtc"
        ),
    output:
        get_pattern_tier(setup, "evt", check_in_cycle=check_in_cycle),
    params:
        timestamp="{timestamp}",
        datatype="{datatype}",
        tier="evt",
        ro_input=lambda _, input: {k: ro(v) for k, v in input.items()},
    log:
        get_pattern_log(setup, f"tier_evt"),
    group:
        "tier-evt"
    resources:
        runtime=300,
        mem_swap=50,
    run:
        shell_string = (
            f"{swenv} python3 -B "
            f"{basedir}/../scripts/build_evt.py "
            f"--configs {ro(configs)} "
            f"--metadata {ro(meta)} "
            "--log {log} "
            "--tier {params.tier} "
            "--datatype {params.datatype} "
            "--timestamp {params.timestamp} "
            "--xtc_file {params.ro_input[xtalk_matrix]} "
            "--par_files {params.ro_input[par_files]} "
            "--hit_file {params.ro_input[hit_file]} "
            "--tcm_file {params.ro_input[tcm_file]} "
            "--dsp_file {params.ro_input[dsp_file]} "
            "--output {output} "
        )
        if input.ann_file is not None:
            shell_string += "--ann_file {params.ro_input[ann_file]} "

        shell(shell_string)


rule build_pet:
    input:
        dsp_file=get_pattern_tier(setup, "psp", check_in_cycle=False),
        hit_file=get_pattern_tier(setup, "pht", check_in_cycle=False),
        tcm_file=get_pattern_tier(setup, "tcm", check_in_cycle=False),
        ann_file=lambda wildcards: (
            None
            if int(wildcards["period"][1:]) > 11
            else get_pattern_tier(setup, "pan", check_in_cycle=False)
        ),
        par_files=lambda wildcards: ParsCatalog.get_par_file(
            setup, wildcards.timestamp, "pht"
        ),
        xtalk_matrix=lambda wildcards: get_input_par_file(
            tier="pet", wildcards=wildcards, name="xtc"
        ),
    output:
        get_pattern_tier(setup, "pet", check_in_cycle=check_in_cycle),
    params:
        timestamp="{timestamp}",
        datatype="{datatype}",
        tier="pet",
        ro_input=lambda _, input: {k: ro(v) for k, v in input.items()},
    log:
        get_pattern_log(setup, f"tier_pet"),
    group:
        "tier-evt"
    resources:
        runtime=300,
        mem_swap=50,
    run:
        shell_string = (
            f"{swenv} python3 -B "
            f"{basedir}/../scripts/build_evt.py "
            f"--configs {ro(configs)} "
            f"--metadata {ro(meta)} "
            "--log {log} "
            "--tier {params.tier} "
            "--datatype {params.datatype} "
            "--timestamp {params.timestamp} "
            "--xtc_file {params.ro_input[xtalk_matrix]} "
            "--par_files {params.ro_input[par_files]} "
            "--hit_file {params.ro_input[hit_file]} "
            "--tcm_file {params.ro_input[tcm_file]} "
            "--dsp_file {params.ro_input[dsp_file]} "
            "--output {output} "
        )
        if input.ann_file is not None:
            shell_string += "--ann_file {params.ro_input[ann_file]} "

        shell(shell_string)


for evt_tier in ("evt", "pet"):

    rule:
        wildcard_constraints:
            timestamp=r"(?!\d{8}T\d{6}Z)",
        input:
            lambda wildcards: sorted(
                get_filelist_full_wildcards(
                    wildcards,
                    setup,
                    get_pattern_tier_raw(setup),
                    tier,
                    ignore_keys_file=os.path.join(configs, "ignore_keys.keylist"),
                )
            ),
        output:
            get_pattern_tier(setup, f"{evt_tier}_concat", check_in_cycle=check_in_cycle),
        params:
            timestamp="all",
            datatype="{datatype}",
            lh5concat_exe=setup["paths"]["install"] + "/bin/lh5concat",
            ro_input=lambda _, input: utils.as_ro(setup, input),
        log:
            get_pattern_log_concat(setup, f"tier_{evt_tier}_concat"),
        group:
            "tier-evt"
        shell:
            "{swenv} {params.lh5concat_exe} --verbose --overwrite "
            "--output {output} "
            "-- {params.ro_input} &> {log}"

    set_last_rule_name(workflow, f"concat_{evt_tier}")
