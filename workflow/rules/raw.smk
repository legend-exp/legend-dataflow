from legenddataflow.patterns import (
    get_pattern_tier_daq_unsorted,
    get_pattern_tier_daq,
    get_pattern_tier,
    get_pattern_log,
    get_pattern_tier_raw_blind,
)
from legenddataflow.paths import (
    config_path,
    chan_map_path,
    metadata_path,
)
from legenddataflow.utils import set_last_rule_name
from legenddataflow.create_pars_keylist import ParsKeyResolve
from legenddataflow.execenv import execenv_pyexe

raw_par_catalog = ParsKeyResolve.get_par_catalog(
    ["-*-*-*-cal"],
    [
        get_pattern_tier_daq_unsorted(config, extension="*"),
        get_pattern_tier_daq(config, extension="*"),
        get_pattern_tier(config, "raw", check_in_cycle=False),
    ],
    {"cal": ["par_raw"]},
)


rule build_raw_orca:
    """
    This rule runs build_raw, it takes in a file.{daq_ext} and outputs a raw file
    """
    input:
        get_pattern_tier_daq(config, extension="orca"),
    params:
        timestamp="{timestamp}",
        datatype="{datatype}",
        ro_input=lambda _, input: ro(input),
        alias_table=lambda wildcards: get_alias(
            channelmap_textdb, wildcards.timestamp, wildcards.datatype, "raw"
        ),
        configs=ro(config_path(config)),
        chan_maps=ro(chan_map_path(config)),
    output:
        get_pattern_tier(config, "raw", check_in_cycle=check_in_cycle),
    log:
        get_pattern_log(config, "tier_raw", time),
    group:
        "tier-raw"
    resources:
        mem_swap=110,
        runtime=300,
    shell:
        execenv_pyexe(config, "build-tier-raw-orca") + "--log {log} "
        "--configs {params.configs} "
        "--chan-maps {params.chan_maps} "
        "--datatype {params.datatype} "
        "--timestamp {params.timestamp} "
        "--alias-table '{params.alias_table}' "
        "{params.ro_input} {output} "


use rule build_raw_orca as build_raw_orca_bz2 with:
    input:
        get_pattern_tier_daq(config, extension="orca.bz2"),


use rule build_raw_orca as build_raw_orca_gz with:
    input:
        get_pattern_tier_daq(config, extension="orca.gz"),


rule build_raw_fcio:
    """
    This rule runs build_raw, it takes in a file.{daq_ext} and outputs a raw file
    """
    input:
        get_pattern_tier_daq(config, extension="fcio"),
    params:
        timestamp="{timestamp}",
        datatype="{datatype}",
        ro_input=lambda _, input: ro(input),
        alias_table=lambda wildcards: get_alias(
            channelmap_textdb, wildcards.timestamp, wildcards.datatype, "raw"
        ),
        configs=ro(config_path(config)),
        chan_maps=ro(chan_map_path(config)),
    output:
        get_pattern_tier(config, "raw", check_in_cycle=check_in_cycle),
    log:
        get_pattern_log(config, "tier_raw", time),
    group:
        "tier-raw"
    resources:
        mem_swap=110,
        runtime=300,
    shell:
        execenv_pyexe(config, "build-tier-raw-fcio") + "--log {log} "
        "--configs {params.configs} "
        "--chan-maps {params.chan_maps} "
        "--datatype {params.datatype} "
        "--timestamp {params.timestamp} "
        "--alias-table '{params.alias_table}' "
        "{params.ro_input} {output}"


rule build_raw_blind:
    """
    This rule runs the data blinding, it takes in the raw file, calibration curve stored in the overrides
    and runs only if the blinding check file is on disk. Output is just the blinded raw file.
    """
    input:
        tier_file=str(get_pattern_tier(config, "raw", check_in_cycle=False)).replace(
            "{datatype}", "phy"
        ),
        blind_file=lambda wildcards: get_blinding_check_file(wildcards, raw_par_catalog),
    params:
        timestamp="{timestamp}",
        datatype="phy",
        ro_input=lambda _, input: {k: ro(v) for k, v in input.items()},
        alias_table=lambda wildcards: get_alias(
            channelmap_textdb, wildcards.timestamp, wildcards.datatype, "raw"
        ),
        configs=ro(config_path(config)),
        chan_maps=ro(chan_map_path(config)),
        metadata=ro(metadata_path(config)),
    output:
        get_pattern_tier_raw_blind(config),
    log:
        str(get_pattern_log(config, "tier_raw_blind", time)).replace(
            "{datatype}", "phy"
        ),
    group:
        "tier-raw"
    resources:
        mem_swap=110,
        runtime=300,
    shell:
        execenv_pyexe(config, "build-tier-raw-blind") + "--log {log} "
        "--configs {params.configs} "
        "--chan-maps {params.chan_maps} "
        "--metadata {params.metadata} "
        "--datatype {params.datatype} "
        "--timestamp {params.timestamp} "
        "--alias-table '{params.alias_table}' "
        "--blind-curve {params.ro_input[blind_file]} "
        "--input {params.ro_input[tier_file]} "
        "--output {output}"
