"""
Snakemake rules for processing psp (partition dsp) tier data.
- combining of all channels into single pars files with associated plot and results files
- running build hit over all channels using par file
"""

from legenddataflow.pars_loading import ParsCatalog
from legenddataflow.create_pars_keylist import ParsKeyResolve
from pathlib import Path
from legenddataflow.patterns import (
    get_pattern_plts,
    get_pattern_tier,
    get_pattern_pars_tmp,
    get_pattern_log,
    get_pattern_pars,
)
from legenddataflow.execenv import execenv_pyexe

psp_par_catalog = ParsKeyResolve.get_par_catalog(
    ["-*-*-*-cal"],
    get_pattern_tier(config, "raw", check_in_cycle=False),
    {"cal": ["par_psp"], "lar": ["par_psp"]},
)


include: "channel_merge.smk"


build_merge_rules("psp", lh5_merge=True, lh5_tier="dsp")


rule build_psp:
    input:
        raw_file=get_pattern_tier(config, "raw", check_in_cycle=False),
        pars_file=ancient(
            lambda wildcards: ParsCatalog.get_par_file(
                psp_par_catalog, config, wildcards.timestamp, "psp"
            )
        ),
    params:
        timestamp="{timestamp}",
        datatype="{datatype}",
        ro_input=lambda _, input: {k: ro(v) for k, v in input.items()},
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
    shell:
        execenv_pyexe(config, "build-tier-dsp") + "--log {log} "
        "--tier psp "
        f"--configs {ro(configs)} "
        "--metadata {meta} "
        "--datatype {params.datatype} "
        "--timestamp {params.timestamp} "
        "--input {params.ro_input[raw_file]} "
        "--output {output.tier_file} "
        "--db-file {output.db_file} "
        "--pars-file {params.ro_input[pars_file]} "
