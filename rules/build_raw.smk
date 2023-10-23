import pathlib, os, json, sys
import scripts as ds
from scripts.util.patterns import *
from datetime import datetime

# Set with `snakemake --configfile=/path/to/your/config.json`
# configfile: "have/to/specify/path/to/your/config.json"

subst_vars_in_snakemake_config(workflow, config)

setup = config["setups"]["l200"]
configs = config_path(setup)
chan_maps = chan_map_path(setup)
swenv = runcmd(setup)

basedir = workflow.basedir


localrules:
    do_nothing,
    gen_filelist,
    autogen_output,


onstart:
    print("Starting workflow")


onsuccess:
    print("Workflow finished, no error")
    shell("rm *.gen || true")
    shell(f"rm {filelist_path(setup)}/* || true")


def get_pattern(tier):
    if tier == "daq":
        return get_pattern_unsorted_data(setup)
    elif tier == "raw":
        return get_pattern_tier_daq(setup)
    else:
        return get_pattern_tier_raw(setup)


checkpoint gen_filelist:
    output:
        os.path.join(filelist_path(setup), "{label}-{tier}.{extension}list"),
    params:
        setup=lambda wildcards: setup,
        search_pattern=lambda wildcards: get_pattern(wildcards.tier),
    script:
        "scripts/create_filelist.py"


def read_filelist(wildcards):
    with checkpoints.gen_filelist.get(
        label=wildcards.label, tier=wildcards.tier, extension="file"
    ).output[0].open() as f:
        files = f.read().splitlines()
        return files


rule autogen_output:
    input:
        filelist=read_filelist,
    output:
        gen_output="{label}-{tier}.gen",
        summary_log=f"{log_path(setup)}/summary-"
        + "{label}-{tier}"
        + f"-{datetime.strftime(datetime.utcnow(), '%Y%m%dT%H%M%SZ')}.log",
        warning_log=f"{log_path(setup)}/warning-"
        + "{label}-{tier}"
        + f"-{datetime.strftime(datetime.utcnow(), '%Y%m%dT%H%M%SZ')}.log",
    params:
        log_path=tmp_log_path(setup),
    shell:
        "{swenv} python3 -B {basedir}/scripts/complete_run.py "
        "--log_path {params.log_path} "
        "--gen_output {output.gen_output} "
        "--summary_log {output.summary_log} "
        "--warning_log {output.warning_log} "
        "--filelist {input.filelist}"


rule sort_data:
    input:
        get_pattern_unsorted_data(setup),
    output:
        get_pattern_tier_daq(setup),
    shell:
        "mv {input} {output}"


rule build_raw:
    input:
        get_pattern_tier_daq(setup),
    params:
        timestamp="{timestamp}",
        datatype="{datatype}",
    output:
        get_pattern_tier_raw(setup),
    log:
        get_pattern_log(setup, "tier_raw"),
    group:
        "tier-raw"
    resources:
        mem_swap=110,
        runtime=300,
    shell:
        "{swenv} python3 -B {basedir}/scripts/build_raw.py "
        "--log {log} "
        "--configs {configs} "
        "--chan_maps {chan_maps} "
        "--datatype {params.datatype} "
        "--timestamp {params.timestamp} "
        "{input} {output}"
