"""
Snakemake rules for processing evt tier.
"""

from scripts.util.pars_loading import ParsCatalog
from scripts.util.patterns import (
    get_pattern_tier,
    get_pattern_log,
    get_pattern_pars,
    get_pattern_log_concat,
)


for tier in ("evt", "pet"):

    rule:
        input:
            dsp_file=(
                get_pattern_tier(setup, "dsp", check_in_cycle=False)
                if tier == "evt"
                else get_pattern_tier(setup, "psp", check_in_cycle=False)
            ),
            hit_file=(
                get_pattern_tier(setup, "hit", check_in_cycle=False)
                if tier == "evt"
                else get_pattern_tier(setup, "pht", check_in_cycle=False)
            ),
            tcm_file=get_pattern_tier(setup, "tcm", check_in_cycle=False),
            xtalk_matrix=lambda wildcards: get_input_par_file(
                tier=tier, wildcards=wildcards, name="xtc"
            ),
            ann_file=branch(
                lambda wildcards: tier if wildcards["period"][1:] <= 11 else "none",
                cases={
                    "evt": get_pattern_tier(setup, "ann", check_in_cycle=False),
                    "pet": get_pattern_tier(setup, "pan", check_in_cycle=False),
                    "none": None,
                },
            ),
            par_files=lambda wildcards: ParsCatalog.get_par_file(
                setup, wildcards.timestamp, "pht"
            ),
        output:
            get_pattern_tier(setup, tier, check_in_cycle=check_in_cycle),
        params:
            timestamp="{timestamp}",
            datatype="{datatype}",
            tier=tier,
            ro_input=lambda _, input: {k: ro(v) for k, v in input.items()},
        log:
            get_pattern_log(setup, f"tier_{tier}"),
        group:
            "tier-evt"
        resources:
            runtime=300,
            mem_swap=50,
        shell:
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
            "--ann_file {params.ro_input[ann_file]} "
            "--output {output} "

    set_last_rule_name(workflow, f"build_{tier}")

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
            get_pattern_tier(setup, f"{tier}_concat", check_in_cycle=check_in_cycle),
        params:
            timestamp="all",
            datatype="{datatype}",
            lh5concat_exe=setup["paths"]["install"] + "/bin/lh5concat",
            ro_input=lambda _, input: utils.as_ro(setup, input),
        log:
            get_pattern_log_concat(setup, f"tier_{tier}_concat"),
        group:
            "tier-evt"
        shell:
            "{swenv} {params.lh5concat_exe} --verbose --overwrite "
            "--output {output} "
            "-- {params.ro_input} &> {log}"

    set_last_rule_name(workflow, f"concat_{tier}")
