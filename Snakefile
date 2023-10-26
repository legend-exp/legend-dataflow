"""
Snakefile for doing the higher stages of data processing (everything beyond build_raw)
This includes:
- building the tcm
- dsp parameter generation
- building dsp
- hit pars generation
- building hit
- building evt
- the same for partition level tiers
"""

import pathlib, os, json, sys
import scripts.util as ds
from scripts.util.pars_loading import pars_catalog
from scripts.util.patterns import *
from scripts.util.utils import (
    subst_vars_in_snakemake_config,
    config_path,
    chan_map_path,
    pars_path,
    filelist_path,
    log_path,
    metadata_path,
    runcmd,
)
from datetime import datetime
from collections import OrderedDict

# Set with `snakemake --configfile=/path/to/your/config.json`
# configfile: "have/to/specify/path/to/your/config.json"

subst_vars_in_snakemake_config(workflow, config)

setup = config["setups"]["l200"]
configs = config_path(setup)
chan_maps = chan_map_path(setup)
meta = metadata_path(setup)
swenv = runcmd(setup)
part = ds.dataset_file(setup, os.path.join(configs, "partitions.json"))
basedir = workflow.basedir


localrules:
    gen_filelist,
    autogen_output,


ds.pars_key_resolve.write_par_catalog(
    ["-*-*-*-cal"],
    os.path.join(pars_path(setup), "pht", "validity.jsonl"),
    get_pattern_tier_raw(setup),
    {"cal": ["par_pht"], "lar": ["par_pht"]},
)


include: "rules/common.smk"
include: "rules/dsp.smk"
include: "rules/hit.smk"
include: "rules/pht.smk"
include: "rules/blinding_check.smk"
include: "rules/blinding_calibration.smk"


onstart:
    print("Starting workflow")
    shell(f"rm {pars_path(setup)}/dsp/validity.jsonl || true")
    shell(f"rm {pars_path(setup)}/hit/validity.jsonl || true")
    shell(f"rm {pars_path(setup)}/pht/validity.jsonl || true")
    shell(f"rm {pars_path(setup)}/raw/validity.jsonl || true")
    ds.pars_key_resolve.write_par_catalog(
        ["-*-*-*-cal"],
        os.path.join(pars_path(setup), "raw", "validity.jsonl"),
        get_pattern_tier_raw(setup),
        {"cal": ["par_raw"]},
    )
    ds.pars_key_resolve.write_par_catalog(
        ["-*-*-*-cal"],
        os.path.join(pars_path(setup), "dsp", "validity.jsonl"),
        get_pattern_tier_raw(setup),
        {"cal": ["par_dsp"], "lar": ["par_dsp"]},
    )
    ds.pars_key_resolve.write_par_catalog(
        ["-*-*-*-cal"],
        os.path.join(pars_path(setup), "hit", "validity.jsonl"),
        get_pattern_tier_raw(setup),
        {"cal": ["par_hit"], "lar": ["par_hit"]},
    )
    ds.pars_key_resolve.write_par_catalog(
        ["-*-*-*-cal"],
        os.path.join(pars_path(setup), "pht", "validity.jsonl"),
        get_pattern_tier_raw(setup),
        {"cal": ["par_pht"], "lar": ["par_pht"]},
    )


onsuccess:
    print("Workflow finished, no error")
    shell("rm *.gen || true")
    shell(f"rm {filelist_path(setup)}/* || true")


# Placeholder, can email or maybe put message in slack
onerror:
    print("An error occurred :( ")


checkpoint gen_filelist:
    output:
        os.path.join(filelist_path(setup), "{label}-{tier}.{extension}list"),
    params:
        setup=lambda wildcards: setup,
        search_pattern=lambda wildcards: get_pattern_tier_raw(setup),
        basedir=basedir,
        configs=configs,
        chan_maps=chan_maps,
        blinding=False,
    script:
        "scripts/create_{wildcards.extension}list.py"


rule gen_fileDB_config:
    output:
        "fdb_config.json",
    script:
        "scripts/gen_fiileDB_config.py"


# Create "{label}-{tier}.gen", based on "{label}.keylist" via
# "{label}-{tier}.filelist". Will implicitly trigger creation of all files
# in "{label}-{tier}.filelist".
# Example: "all[-{detector}[-{measurement}[-{run}[-{timestamp}]]]]-{tier}.gen":
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
        tmp_par_path=os.path.join(tmp_par_path(setup), "*_db.json"),
        valid_keys_path=os.path.join(pars_path(setup), "valid_keys"),
        filedb_path=os.path.join(pars_path(setup), "filedb"),
        setup=lambda wildcards: setup,
        basedir=basedir,
    script:
        "scripts/complete_run.py"
