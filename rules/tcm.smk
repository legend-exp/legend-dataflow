"""
Snakemake file containing the rules for generating the tcm
"""

# This rule builds the tcm files each raw file
rule build_tier_tcm:
    input:
        get_pattern_tier_raw(setup),
    output:
        get_pattern_tier_tcm(setup),
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