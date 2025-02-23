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

# FIXME: this seems quite stupid
dtypes_w_pars = ("lac", "cal", "acs", "anc", "anp", "ant", "aph", "ath", "tst")

dsp_par_catalog = ParsKeyResolve.get_par_catalog(
    [f"-*-*-*-{dt}" for dt in dtypes_w_pars],
    patt.get_pattern_tier(config, "raw", check_in_cycle=False),
    {dt: ["par_dsp"] for dt in dtypes_w_pars},
)

build_merge_rules("dsp", lh5_merge=True)


rule build_dsp:
    input:
        raw_file=patt.get_pattern_tier(config, "raw", check_in_cycle=False),
        pars_file=ancient(
            lambda wildcards: dsp_par_catalog.get_par_file(
                config, wildcards.timestamp, "dsp"
            )
        ),
    params:
        timestamp="{timestamp}",
        datatype="{datatype}",
        ro_input=lambda _, input: {k: ro(v) for k, v in input.items()},
    output:
        tier_file=patt.get_pattern_tier(config, "dsp", check_in_cycle=check_in_cycle),
        db_file=patt.get_pattern_pars_tmp(config, "dsp_db"),
    log:
        patt.get_pattern_log(config, "tier_dsp", time),
    group:
        "tier-dsp"
    resources:
        runtime=300,
        mem_swap=lambda wildcards: 35 if wildcards.datatype == "cal" else 25,
    shell:
        execenv_pyexe(config, "build-tier-dsp") + "--log {log} "
        "--tier dsp "
        f"--configs {ro(configs)} "
        "--metadata {meta} "
        "--datatype {params.datatype} "
        "--timestamp {params.timestamp} "
        "--input {params.ro_input[raw_file]} "
        "--output {output.tier_file} "
        "--db-file {output.db_file} "
        "--pars-file {params.ro_input[pars_file]} "
