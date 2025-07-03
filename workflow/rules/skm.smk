"""
Snakemake rules for processing skm tier.
"""

from legenddataflow.methods.patterns import (
    get_pattern_tier,
    get_pattern_log,
    get_pattern_pars,
    get_pattern_log_concat,
)
from legenddataflow.methods.paths import config_path
from legenddataflowscripts.workflow import execenv_pyexe


rule build_skm:
    input:
        get_pattern_tier(config, "pet_concat", check_in_cycle=False),
    output:
        get_pattern_tier(config, "skm", check_in_cycle=check_in_cycle),
    params:
        timestamp="20230410T000000Z",
        datatype="phy",
        ro_input=lambda _, input: ro(input),
        configs=ro(config_path(config)),
    log:
        get_pattern_log_concat(config, "tier_skm", time),
    group:
        "tier-skm"
    resources:
        runtime=300,
    shell:
        execenv_pyexe(config, "build-tier-skm") + "--configs {params.configs} "
        "--timestamp {params.timestamp} "
        "--log {log} "
        "--datatype {params.datatype} "
        "--evt-file {params.ro_input} "
        "--output {output} "
