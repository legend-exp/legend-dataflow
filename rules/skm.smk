"""
Snakemake rules for processing skm tier.
"""

from scripts.util.patterns import (
    get_pattern_tier,
    get_pattern_log,
    get_pattern_pars,
    get_pattern_log_concat,
)


rule build_skm:
    input:
        evt_file = get_pattern_tier(setup, "pet_concat", check_in_cycle=False),
    output:
        skm_file=get_pattern_tier(setup, "skm", check_in_cycle=check_in_cycle),
    params:
        timestamp="all",
        datatype="phy",
    log:
        get_pattern_log_concat(setup, "tier_skm"),
    group:
        "tier-skm"
    resources:
        runtime=300,
    shell:
        "{swenv} python3 -B "
        f"{basedir}/../scripts/build_skm.py "
        "--configs {configs} "
        "--metadata {meta} "
        "--log {log} "
        "--datatype {params.datatype} "
        "--evt_file {input.evt_file} "
        "--output {output.skm_file} "
