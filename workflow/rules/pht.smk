"""
Snakemake rules for processing pht (partition hit) tier data. This is done in 4 steps:
- extraction of calibration curves(s) for each run for each channel from cal data
- extraction of psd calibration parameters and partition level energy fitting for each channel over whole partition from cal data
- combining of all channels into single pars files with associated plot and results files
- running build hit over all channels using par file
"""

from pathlib import Path
from legenddataflow.methods import ParsCatalog, ParsKeyResolve
from legenddataflow.methods.paths import filelist_path, config_path, metadata_path
from legenddataflow.methods.patterns import (
    get_pattern_tier,
    get_pattern_pars_tmp,
    get_pattern_log,
)
from legenddataflowscripts.workflow import execenv_pyexe, set_last_rule_name

pht_par_catalog = ParsKeyResolve.get_par_catalog(
    ["-*-*-*-cal"],
    get_pattern_tier(config, "raw", check_in_cycle=False),
    {"cal": ["par_pht"], "lar": ["par_pht"]},
    run_overwrite_validity,
    ignore_keys_file=Path(det_status) / "ignored_daq_cycles.yaml",
)

intier = config.get("pht_intier", "psp")

build_merge_rules("pht", lh5_merge=False)


rule build_pht:
    input:
        dsp_file=get_pattern_tier(config, intier, check_in_cycle=False),
        pars_file=lambda wildcards: pht_par_catalog.get_par_file(
            config, wildcards.timestamp, "pht"
        ),
    output:
        tier_file=get_pattern_tier(config, "pht", check_in_cycle=check_in_cycle),
    params:
        timestamp="{timestamp}",
        datatype="{datatype}",
        tier="pht",
        ro_input=lambda _, input: {k: ro(v) for k, v in input.items()},
        table_map=lambda wildcards: get_table_mapping(
            channelmap_textdb, wildcards.timestamp, wildcards.datatype, "dsp"
        ),
        alias_table=lambda wildcards: get_alias(
            channelmap_textdb, wildcards.timestamp, wildcards.datatype, "hit"
        ),
        configs=ro(config_path(config)),
    log:
        get_pattern_log(config, "tier_pht", time),
    group:
        "tier-pht"
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
