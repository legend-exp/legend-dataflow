"""
Snakemake rules for processing ann tier. This is done only for the coax detectors
to apply the ann and risetime cuts for psd.

"""

from legenddataflow.patterns import (
    get_pattern_tier,
    get_pattern_log,
    get_pattern_pars,
)


rule build_ann:
    input:
        dsp_file=get_pattern_tier(setup, "dsp", check_in_cycle=False),
        pars_file=lambda wildcards: get_input_par_file(wildcards, "ann", "cuts"),
    params:
        timestamp="{timestamp}",
        datatype="{datatype}",
    output:
        tier_file=get_pattern_tier(setup, "ann", check_in_cycle=check_in_cycle),
        db_file=get_pattern_pars_tmp(setup, "ann_db"),
    log:
        get_pattern_log(setup, "tier_ann", time),
    group:
        "tier-ann"
    resources:
        runtime=300,
        mem_swap=lambda wildcards: 25 if wildcards.datatype == "cal" else 15,
    shell:
        "{swenv} python3 -B "
        f"{workflow.source_path('../scripts/build_dsp.py')} "
        "--log {log} "
        "--configs {configs} "
        "--metadata {meta} "
        f"--tier ann "
        "--datatype {params.datatype} "
        "--timestamp {params.timestamp} "
        "--input {input.dsp_file} "
        "--output {output.tier_file} "
        "--db_file {output.db_file} "
        "--pars_file {input.pars_file} "


rule build_pan:
    input:
        dsp_file=get_pattern_tier(setup, "psp", check_in_cycle=False),
        pars_file=lambda wildcards: get_input_par_file(wildcards, "ann", "cuts"),
    params:
        timestamp="{timestamp}",
        datatype="{datatype}",
    output:
        tier_file=get_pattern_tier(setup, "pan", check_in_cycle=check_in_cycle),
        db_file=get_pattern_pars_tmp(setup, "pan_db"),
    log:
        get_pattern_log(setup, "tier_pan", time),
    group:
        "tier-ann"
    resources:
        runtime=300,
        mem_swap=lambda wildcards: 25 if wildcards.datatype == "cal" else 15,
    shell:
        "{swenv} python3 -B "
        f"{workflow.source_path('../scripts/build_dsp.py')} "
        "--log {log} "
        "--configs {configs} "
        "--metadata {meta} "
        f"--tier pan "
        "--datatype {params.datatype} "
        "--timestamp {params.timestamp} "
        "--input {input.dsp_file} "
        "--output {output.tier_file} "
        "--db_file {output.db_file} "
        "--pars_file {input.pars_file} "
