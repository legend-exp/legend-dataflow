"""
Snakemake rules for processing ann tier. This is done only for the coax detectors
to apply the ann and risetime cuts for psd.

"""

from legenddataflow.patterns import (
    get_pattern_tier,
    get_pattern_log,
    get_pattern_pars,
)
from legenddataflow.paths import config_path
from legenddataflow.execenv import execenv_pyexe


rule build_ann:
    input:
        dsp_file=get_pattern_tier(config, "dsp", check_in_cycle=False),
        pars_file=lambda wildcards: get_input_par_file(
            config=config, wildcards=wildcards, tier="ann", name="cuts"
        ),
    params:
        timestamp="{timestamp}",
        datatype="{datatype}",
        table_map=lambda wildcards: get_table_mapping(
            channelmap_textdb, wildcards.timestamp, wildcards.datatype, "dsp"
        ),
        alias_table=lambda wildcards: get_alias(
            channelmap_textdb, wildcards.timestamp, wildcards.datatype, "dsp"
        ),
        configs=ro(config_path(config)),
        meta=ro(metadata_path(config)),
    output:
        tier_file=get_pattern_tier(config, "ann", check_in_cycle=check_in_cycle),
    log:
        get_pattern_log(config, "tier_ann", time),
    group:
        "tier-ann"
    resources:
        runtime=300,
        mem_swap=lambda wildcards: 25 if wildcards.datatype == "cal" else 15,
    threads: 5 if config.get("multiprocess", False) else 1
    shell:
        execenv_pyexe(config, "build-tier-dsp") + "--log {log} "
        "--configs {params.configs} "
        "--table-map '{params.table_map}' "
        "--tier ann "
        "--datatype {params.datatype} "
        "--timestamp {params.timestamp} "
        "--alias-table '{params.alias_table}' "
        "--input {input.dsp_file} "
        "--output {output.tier_file} "
        "--pars-file {input.pars_file} "
        "--n-processes {threads} "


rule build_pan:
    input:
        dsp_file=get_pattern_tier(config, "psp", check_in_cycle=False),
        pars_file=lambda wildcards: get_input_par_file(
            config=config, wildcards=wildcards, tier="ann", name="cuts"
        ),
    params:
        timestamp="{timestamp}",
        datatype="{datatype}",
        table_map=lambda wildcards: get_table_mapping(
            channelmap_textdb, wildcards.timestamp, wildcards.datatype, "dsp"
        ),
        alias_table=lambda wildcards: get_alias(
            channelmap_textdb, wildcards.timestamp, wildcards.datatype, "dsp"
        ),
        configs=ro(config_path(config)),
    output:
        tier_file=get_pattern_tier(config, "pan", check_in_cycle=check_in_cycle),
    log:
        get_pattern_log(config, "tier_pan", time),
    group:
        "tier-ann"
    resources:
        runtime=300,
        mem_swap=lambda wildcards: 25 if wildcards.datatype == "cal" else 15,
    threads: 5 if config.get("multiprocess", False) else 1
    shell:
        execenv_pyexe(config, "build-tier-dsp") + "--log {log} "
        "--configs {params.configs} "
        "--table-map '{params.table_map}' "
        "--tier pan "
        "--datatype {params.datatype} "
        "--timestamp {params.timestamp} "
        "--alias-table '{params.alias_table}' "
        "--input {input.dsp_file} "
        "--output {output.tier_file} "
        "--pars-file {input.pars_file} "
        "--n-processes {threads} "
