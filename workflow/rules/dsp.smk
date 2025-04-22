"""
Snakemake rules for processing dsp tier.
- combining of all channels into single pars files with associated plot and results files
- running dsp over all channels using par file
"""

from legenddataflow.pars_loading import ParsCatalog
from legenddataflow.create_pars_keylist import ParsKeyResolve
from pathlib import Path
from legenddataflow import patterns as patt
from legenddataflow.execenv import execenv_pyexe

build_merge_rules("dsp", lh5_merge=True)

dsp_par_catalog = ParsKeyResolve.get_par_catalog(
    ["-*-*-*-cal"],
    get_pattern_tier(config, "raw", check_in_cycle=False),
    {"cal": ["par_dsp"], "lar": ["par_dsp"]},
)


rule build_dsp:
    input:
        raw_file=patt.get_pattern_tier(config, "raw", check_in_cycle=False),
        pars_files=ancient(lambda wildcards: _make_input_pars_file(wildcards)),
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
        tier_file=patt.get_pattern_tier(config, "dsp", check_in_cycle=check_in_cycle),
    log:
        patt.get_pattern_log(config, "tier_dsp", time),
    group:
        "tier-dsp"
    resources:
        runtime=300,
        mem_swap=lambda wildcards: 90 if wildcards.datatype == "cal" else 60,
    threads: get_threads
    shell:
        execenv_pyexe(config, "build-tier-dsp") + "--log {log} "
        "--tier dsp "
        "--configs {params.configs} "
        "--table-map '{params.table_map}' "
        "--datatype {params.datatype} "
        "--timestamp {params.timestamp} "
        "--alias-table '{params.alias_table}' "
        "--input {params.ro_input[raw_file]} "
        "--output {output.tier_file} "
        "--pars-file {params.ro_input[pars_files]} "
        "--n-processes {threads} "
