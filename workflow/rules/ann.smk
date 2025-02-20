"""
Snakemake rules for processing ann tier. This is done only for the coax detectors
to apply the ann and risetime cuts for psd.

"""

from legenddataflow.patterns import (
    get_pattern_tier,
    get_pattern_log,
    get_pattern_pars,
)
from legenddataflow.execenv import execenv_pyexe


rule build_ann:
    input:
        dsp_file=get_pattern_tier(config, "dsp", check_in_cycle=False),
        pars_file=lambda wildcards: get_input_par_file(
            setup=config, wildcards=wildcards, tier="ann", name="cuts"
        ),
    params:
        timestamp="{timestamp}",
        datatype="{datatype}",
    output:
        tier_file=get_pattern_tier(config, "ann", check_in_cycle=check_in_cycle),
        db_file=get_pattern_pars_tmp(config, "ann_db"),
    log:
        get_pattern_log(config, "tier_ann", time),
    group:
        "tier-ann"
    resources:
        runtime=300,
        mem_swap=lambda wildcards: 25 if wildcards.datatype == "cal" else 15,
    shell:
        execenv_pyexe(config, "build-tier-dsp") + "--log {log} "
        "--configs {configs} "
        "--metadata {meta} "
        "--tier ann "
        "--datatype {params.datatype} "
        "--timestamp {params.timestamp} "
        "--input {input.dsp_file} "
        "--output {output.tier_file} "
        "--db-file {output.db_file} "
        "--pars-file {input.pars_file} "


rule build_pan:
    input:
        dsp_file=get_pattern_tier(config, "psp", check_in_cycle=False),
        pars_file=lambda wildcards: get_input_par_file(
            setup=config, wildcards=wildcards, tier="ann", name="cuts"
        ),
    params:
        timestamp="{timestamp}",
        datatype="{datatype}",
    output:
        tier_file=get_pattern_tier(config, "pan", check_in_cycle=check_in_cycle),
        db_file=get_pattern_pars_tmp(config, "pan_db"),
    log:
        get_pattern_log(config, "tier_pan", time),
    group:
        "tier-ann"
    resources:
        runtime=300,
        mem_swap=lambda wildcards: 25 if wildcards.datatype == "cal" else 15,
    shell:
        execenv_pyexe(config, "build-tier-dsp") + "--log {log} "
        "--configs {configs} "
        "--metadata {meta} "
        "--tier pan "
        "--datatype {params.datatype} "
        "--timestamp {params.timestamp} "
        "--input {input.dsp_file} "
        "--output {output.tier_file} "
        "--db-file {output.db_file} "
        "--pars-file {input.pars_file} "
