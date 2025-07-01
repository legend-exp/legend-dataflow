"""
Snakemake rules for processing hit tier. This is done in 4 steps:
- extraction of calibration curves(s) for each channel from cal data
- extraction of psd calibration parameters for each channel from cal data
- combining of all channels into single pars files with associated plot and results files
- running build hit over all channels using par file
"""

from legenddataflow.methods.pars_loading import ParsCatalog, ParsKeyResolve
from pathlib import Path
from legenddataflow.methods.patterns import (
    get_pattern_tier,
    get_pattern_log,
    get_pattern_pars_tmp,
)
from legenddataflow.methods.paths import config_path
from legenddataflowscripts.workflow import execenv_pyexe

hit_par_catalog = ParsKeyResolve.get_par_catalog(
    ["-*-*-*-cal"],
    get_pattern_tier(config, "raw", check_in_cycle=False),
    {"cal": ["par_hit"], "lar": ["par_hit"]},
)

build_merge_rules("hit", lh5_merge=False)


rule build_hit:
    input:
        dsp_file=get_pattern_tier(config, "dsp", check_in_cycle=False),
        pars_file=lambda wildcards: hit_par_catalog.get_par_file(
            config, wildcards.timestamp, "hit"
        ),
    output:
        tier_file=get_pattern_tier(config, "hit", check_in_cycle=check_in_cycle),
    params:
        timestamp="{timestamp}",
        datatype="{datatype}",
        tier="hit",
        ro_input=lambda _, input: {k: ro(v) for k, v in input.items()},
        table_map=lambda wildcards: get_table_mapping(
            channelmap_textdb, wildcards.timestamp, wildcards.datatype, "dsp"
        ),
        alias_table=lambda wildcards: get_alias(
            channelmap_textdb, wildcards.timestamp, wildcards.datatype, "hit"
        ),
        configs=ro(config_path(config)),
    log:
        get_pattern_log(config, "tier_hit", time),
    group:
        "tier-hit"
    resources:
        runtime=300,
    shell:
        execenv_pyexe(config, "build-tier-hit") + "--configs {params.configs} "
        "--table-map '{params.table_map}' "
        "--tier {params.tier} "
        "--datatype {params.datatype} "
        "--timestamp {params.timestamp} "
        "--alias-table '{params.alias_table}' "
        "--pars-file {params.ro_input[pars_file]} "
        "--output {output.tier_file} "
        "--input {params.ro_input[dsp_file]} "
        "--log {log}"
