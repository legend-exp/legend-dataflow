"""
Snakemake rules for processing hit tier. This is done in 4 steps:
- extraction of calibration curves(s) for each channel from cal data
- extraction of psd calibration parameters for each channel from cal data
- combining of all channels into single pars files with associated plot and results files
- running build hit over all channels using par file
"""

from legenddataflow.create_pars_keylist import ParsKeyResolve
from legenddataflow.pars_loading import ParsCatalog
from pathlib import Path
from legenddataflow.patterns import (
    get_pattern_tier,
    get_pattern_log,
    get_pattern_pars_tmp,
)
from legenddataflow.execenv import execenv_smk_py_script

hit_par_catalog = ParsKeyResolve.get_par_catalog(
    ["-*-*-*-cal"],
    get_pattern_tier(config, "raw", check_in_cycle=False),
    {"cal": ["par_hit"], "lar": ["par_hit"]},
)


include: "channel_merge.smk"


build_merge_rules("hit", lh5_merge=False)


rule build_hit:
    input:
        dsp_file=get_pattern_tier(config, "dsp", check_in_cycle=False),
        pars_file=lambda wildcards: ParsCatalog.get_par_file(
            config, wildcards.timestamp, "hit"
        ),
    output:
        tier_file=get_pattern_tier(config, "hit", check_in_cycle=check_in_cycle),
        db_file=get_pattern_pars_tmp(config, "hit_db"),
    params:
        timestamp="{timestamp}",
        datatype="{datatype}",
        tier="hit",
        ro_input=lambda _, input: {k: ro(v) for k, v in input.items()},
    log:
        get_pattern_log(config, "tier_hit", time),
    group:
        "tier-hit"
    resources:
        runtime=300,
    shell:
        f'{execenv_smk_py_script(config, "build-tier-hit")}'
        f"--configs {ro(configs)} "
        "--metadata {meta} "
        "--log {log} "
        "--tier {params.tier} "
        "--datatype {params.datatype} "
        "--timestamp {params.timestamp} "
        "--pars-file {params.ro_input[pars_file]} "
        "--output {output.tier_file} "
        "--input {params.ro_input[dsp_file]} "
        "--db-file {output.db_file}"
