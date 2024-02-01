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


rule build_evt:
    input:
        dsp_file=get_pattern_tier_dsp(setup),
        hit_file=get_pattern_tier_hit(setup),
        tcm_file=get_pattern_tier_tcm(setup),
    output:
        evt_file=get_pattern_tier(setup, "evt", check_in_cycle=check_in_cycle),
    params:
        timestamp="{timestamp}",
        datatype="{datatype}",
        tier="evt",
    log:
        get_pattern_log(setup, "tier_evt"),
    group:
        "tier-evt"
    resources:
        runtime=300,
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


rule build_pet:
    input:
        dsp_file=get_pattern_tier_dsp(setup),
        hit_file=get_pattern_tier_pht(setup),
        tcm_file=get_pattern_tier_tcm(setup),
    output:
        evt_file=get_pattern_tier(setup, "pet", check_in_cycle=check_in_cycle),
    params:
        timestamp="{timestamp}",
        datatype="{datatype}",
        tier="pet",
    log:
        get_pattern_log(setup, "tier_pet"),
    group:
        "tier-evt"
    resources:
        runtime=300,
    shell:
        "{swenv} python3 -B "
        f"{workflow.source_path('../scripts/build_evt.py')} "
        "--configs {configs} "
        "--log {log} "
        "--tier {params.tier} "
        "--datatype {params.datatype} "
        "--timestamp {params.timestamp} "
        "--metadata {meta} "
        "--hit_file {input.hit_file} "
        "--tcm_file {input.tcm_file} "
        "--dsp_file {input.dsp_file} "
        "--output {output.evt_file} "
