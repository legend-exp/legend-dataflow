"""
Snakemake file containing the rules for generating the tcm
"""

from scripts.util.patterns import (
    get_pattern_tier_raw,
    get_pattern_tier,
    get_pattern_log,
)


# This rule builds the tcm files each raw file
rule build_tier_tcm:
    input:
        get_pattern_tier_raw(setup),
    output:
        get_pattern_tier(setup, "tcm", check_in_cycle=check_in_cycle),
    log:
        get_pattern_log(setup, "tier_tcm"),
    group:
        "tier-tcm"
    resources:
        runtime=300,
    shell:
        "{swenv} python3 -B "
        f"{workflow.source_path('../scripts/build_tcm.py')} "
        "--log {log} "
        "--configs {configs} "
        "{input} "
        "{output}"
