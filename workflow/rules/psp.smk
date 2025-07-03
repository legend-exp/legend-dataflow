"""
Snakemake rules for processing psp (partition dsp) tier data.
- combining of all channels into single pars files with associated plot and results files
- running build hit over all channels using par file
"""

from pathlib import Path
from legenddataflow.methods import ParsKeyResolve, ParsCatalog
from legenddataflow.methods.patterns import (
    get_pattern_plts,
    get_pattern_tier,
    get_pattern_pars_tmp,
    get_pattern_log,
    get_pattern_pars,
)
from legenddataflow.methods.paths import config_path
from legenddataflowscripts.workflow import execenv_pyexe

psp_par_catalog = ParsKeyResolve.get_par_catalog(
    ["-*-*-*-cal"],
    get_pattern_tier(config, "raw", check_in_cycle=False),
    {"cal": ["par_psp"], "lar": ["par_psp"]},
)

build_merge_rules("psp", lh5_merge=True, lh5_tier="dsp")


rule build_psp:
    input:
        raw_file=get_pattern_tier(config, "raw", check_in_cycle=False),
        pars_file=ancient(
            lambda wildcards: psp_par_catalog.get_par_file(
                config, wildcards.timestamp, "psp"
            )
        ),
    params:
        timestamp="{timestamp}",
        datatype="{datatype}",
        ro_input=lambda _, input: {k: ro(v) for k, v in input.items()},
        table_map=lambda wildcards: get_table_mapping(
            channelmap_textdb, wildcards.timestamp, wildcards.datatype, "raw"
        ),
        alias_table=lambda wildcards: get_alias(
            channelmap_textdb, wildcards.timestamp, wildcards.datatype, "dsp"
        ),
        configs=ro(config_path(config)),
    output:
        tier_file=get_pattern_tier(config, "psp", check_in_cycle=check_in_cycle),
        db_file=get_pattern_pars_tmp(config, "psp_db"),
    log:
        get_pattern_log(config, "tier_psp", time),
    group:
        "tier-dsp"
    resources:
        runtime=300,
        mem_swap=lambda wildcards: 35 if wildcards.datatype == "cal" else 25,
    threads: get_threads
    shell:
        execenv_pyexe(config, "build-tier-dsp") + "--log {log} "
        "--tier psp "
        "--configs {params.configs} "
        "--table-map '{params.table_map}' "
        "--alias-table '{params.alias_table}' "
        "--metadata {meta} "
        "--datatype {params.datatype} "
        "--timestamp {params.timestamp} "
        "--input {params.ro_input[raw_file]} "
        "--output {output.tier_file} "
        "--db-file {output.db_file} "
        "--pars-file {params.ro_input[pars_file]} "
        "--n-processes {threads} "
