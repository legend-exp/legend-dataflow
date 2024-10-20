"""
Snakemake rules for processing ann tier. This is done only for the coax detectors
to apply the ann and risetime cuts for psd.

"""

from scripts.util.pars_loading import pars_catalog
from scripts.util.utils import par_dsp_path
from scripts.util.patterns import (
    get_pattern_tier_dsp,
    get_pattern_tier_psp,
    get_pattern_tier_ann,
    get_pattern_tier,
    get_pattern_log,
    get_pattern_pars,
    get_pattern_pars_overwrite,
)

for tier in ["ann", "pan"]:

    rule:
        input:
            dsp_file=get_pattern_tier_dsp(setup) if tier == "ann" else get_pattern_tier_psp(setup),
            pars_file=lambda wildcards: get_svm_file(wildcards, "ann", "cuts"),
        params:
            timestamp="{timestamp}",
            datatype="{datatype}",
        output:
            tier_file=get_pattern_tier(setup, tier, check_in_cycle=check_in_cycle),
            db_file=get_pattern_pars_tmp(setup, f"{tier}_db"),
        log:
            get_pattern_log(setup, f"tier_{tier}"),
        group:
            "tier-ann"
        resources:
            runtime=300,
            mem_swap=lambda wildcards: 25 if wildcards.datatype == "cal" else 15,
        shell:
            "{swenv} python3 -B "
            f"{workflow.source_path('../scripts/build_ann.py')} "
            "--log {log} "
            "--configs {configs} "
            "--datatype {params.datatype} "
            "--timestamp {params.timestamp} "
            "--input {input.dsp_file} "
            "--output {output.tier_file} "
            "--db_file {output.db_file} "
            "--pars_file {input.pars_file} "
    
    set_last_rule_name(workflow, f"build_{tier}")