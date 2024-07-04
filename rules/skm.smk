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
        get_pattern_tier(setup, "pet_concat", check_in_cycle=False),
    output:
        get_pattern_tier(setup, "skm", check_in_cycle=check_in_cycle),
    params:
        timestamp="20230410T000000Z",
        datatype="phy",
        ro_input=lambda _, input: ro(input),
    log:
        get_pattern_log_concat(setup, "tier_skm"),
    group:
        "tier-skm"
    resources:
        runtime=300,
    shell:
        "{swenv} python3 -B "
        "{basedir}/../scripts/build_skm.py "
        f"--configs {ro(configs)} "
        "--timestamp {params.timestamp} "
        "--log {log} "
        "--datatype {params.datatype} "
        "--evt_file {params.ro_input} "
        "--output {output} "
