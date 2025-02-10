"""
Snakemake rules for processing pht (partition hit) tier data. This is done in 4 steps:
- extraction of calibration curves(s) for each run for each channel from cal data
- extraction of psd calibration parameters and partition level energy fitting for each channel over whole partition from cal data
- combining of all channels into single pars files with associated plot and results files
- running build hit over all channels using par file
"""

from legenddataflow.create_pars_keylist import ParsKeyResolve
from legenddataflow.pars_loading import ParsCatalog
from pathlib import Path
from legenddataflow.utils import filelist_path, set_last_rule_name
from legenddataflow.patterns import (
    get_pattern_tier,
    get_pattern_pars_tmp,
    get_pattern_log,
)
from legenddataflow.execenv import execenv_pyexe

pht_par_catalog = ParsKeyResolve.get_par_catalog(
    ["-*-*-*-cal"],
    get_pattern_tier(config, "raw", check_in_cycle=False),
    {"cal": ["par_pht"], "lar": ["par_pht"]},
)

intier = "psp"


include: "channel_merge.smk"


build_merge_rules("pht", lh5_merge=False)


rule build_pht:
    input:
        dsp_file=get_pattern_tier(config, intier, check_in_cycle=False),
        pars_file=lambda wildcards: pht_par_catalog.get_par_file(
            config, wildcards.timestamp, "pht"
        ),
    output:
        tier_file=get_pattern_tier(config, "pht", check_in_cycle=check_in_cycle),
        db_file=get_pattern_pars_tmp(config, "pht_db"),
    params:
        timestamp="{timestamp}",
        datatype="{datatype}",
        tier="pht",
        ro_input=lambda _, input: {k: ro(v) for k, v in input.items()},
    log:
        get_pattern_log(config, "tier_pht", time),
    group:
        "tier-pht"
    resources:
        runtime=300,
    shell:
        execenv_pyexe(config, "build-tier-hit") + f"--configs {ro(configs)} "
        "--metadata {meta} "
        "--log {log} "
        "--tier {params.tier} "
        "--datatype {params.datatype} "
        "--timestamp {params.timestamp} "
        "--pars-file {params.ro_input[pars_file]} "
        "--output {output.tier_file} "
        "--input {params.ro_input[dsp_file]} "
        "--db-file {output.db_file}"
