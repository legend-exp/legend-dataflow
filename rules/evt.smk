"""
Snakemake rules for processing evt tier.
"""

from scripts.util.patterns import (
    get_pattern_tier_hit,
    get_pattern_tier_dsp,
    get_pattern_tier_tcm,
    get_pattern_tier_pht,
    get_pattern_tier_psp,
    get_pattern_tier,
    get_pattern_log,
    get_pattern_pars,
)


for tier in ("evt", "pet"):

    rule:
        input:
            dsp_file=(
                get_pattern_tier_dsp(setup)
                if tier == "evt"
                else get_pattern_tier_psp(setup)
            ),
            hit_file=(
                get_pattern_tier_hit(setup)
                if tier == "evt"
                else get_pattern_tier_pht(setup)
            ),
            tcm_file=get_pattern_tier_tcm(setup),
        output:
            evt_file=get_pattern_tier(setup, tier, check_in_cycle=check_in_cycle),
        params:
            timestamp="{timestamp}",
            datatype="{datatype}",
            tier=tier,
        log:
            get_pattern_log(setup, f"tier_{tier}"),
        group:
            "tier-evt"
        resources:
            runtime=300,
            mem_swap=70,
        shell:
            "{swenv} python3 -B "
            f"{workflow.source_path('../scripts/build_evt.py')} "
            "--configs {configs} "
            "--metadata {meta} "
            "--log {log} "
            "--tier {params.tier} "
            "--datatype {params.datatype} "
            "--timestamp {params.timestamp} "
            "--hit_file {input.hit_file} "
            "--tcm_file {input.tcm_file} "
            "--dsp_file {input.dsp_file} "
            "--output {output.evt_file} "

    set_last_rule_name(workflow, f"build_{tier}")

    rule:
        wildcard_constraints:
            timestamp="(?!\d{8}T\d{6}Z)",
        input:
            lambda wildcards: sorted(read_filelist_phy(wildcards, tier)),
        output:
            get_pattern_tier(setup, tier, check_in_cycle=check_in_cycle),
        params:
            timestamp="all",
            datatype="{datatype}",
        log:
            get_pattern_log(setup, "tier_skm"),
        group:
            "tier-evt"
        shell:
            "{swenv} lh5concat --verbose --overwrite "
            "--output {output} "
            "-- {input} &> {log}"

    set_last_rule_name(workflow, f"concat_{tier}")
