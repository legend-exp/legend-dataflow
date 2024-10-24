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

import pathlib
import os
import json
import sys
import glob
from datetime import datetime
from collections import OrderedDict
import logging

import scripts.util as ds
from scripts.util.pars_loading import pars_catalog
from scripts.util.patterns import get_pattern_tier_raw
from scripts.util.utils import (
    subst_vars_in_snakemake_config,
    runcmd,
    config_path,
    chan_map_path,
    filelist_path,
    metadata_path,
    tmp_log_path,
    pars_path,
)

# Set with `snakemake --configfile=/path/to/your/config.json`
# configfile: "have/to/specify/path/to/your/config.json"

subst_vars_in_snakemake_config(workflow, config)

check_in_cycle = True
setup = config["setups"]["l200"]
configs = config_path(setup)
chan_maps = chan_map_path(setup)
meta = metadata_path(setup)
swenv = runcmd(setup)
part = ds.dataset_file(setup, os.path.join(configs, "partitions.json"))
basedir = workflow.basedir


wildcard_constraints:
    experiment=r"\w+",
    period=r"p\d{2}",
    run=r"r\d{3}",
    datatype=r"\w{3}",
    timestamp=r"\d{8}T\d{6}Z",


include: "rules/filelist_gen.smk"
include: "rules/chanlist_gen.smk"
include: "rules/common.smk"
include: "rules/main.smk"
include: "rules/tcm.smk"
include: "rules/dsp.smk"
include: "rules/psp.smk"
include: "rules/hit.smk"
include: "rules/pht.smk"
include: "rules/pht_fast.smk"
include: "rules/evt.smk"
include: "rules/skm.smk"
include: "rules/blinding_calibration.smk"
include: "rules/qc_phy.smk"


# Log parameter catalogs in validity.jsonl files
hit_par_cat_file = os.path.join(pars_path(setup), "hit", "validity.jsonl")
if os.path.isfile(hit_par_cat_file):
    os.remove(os.path.join(pars_path(setup), "hit", "validity.jsonl"))
pathlib.Path(os.path.dirname(hit_par_cat_file)).mkdir(parents=True, exist_ok=True)
ds.pars_key_resolve.write_to_jsonl(hit_par_catalog, hit_par_cat_file)

pht_par_cat_file = os.path.join(pars_path(setup), "pht", "validity.jsonl")
if os.path.isfile(pht_par_cat_file):
    os.remove(os.path.join(pars_path(setup), "pht", "validity.jsonl"))
pathlib.Path(os.path.dirname(pht_par_cat_file)).mkdir(parents=True, exist_ok=True)
ds.pars_key_resolve.write_to_jsonl(pht_par_catalog, pht_par_cat_file)

dsp_par_cat_file = os.path.join(pars_path(setup), "dsp", "validity.jsonl")
if os.path.isfile(dsp_par_cat_file):
    os.remove(dsp_par_cat_file)
pathlib.Path(os.path.dirname(dsp_par_cat_file)).mkdir(parents=True, exist_ok=True)
ds.pars_key_resolve.write_to_jsonl(dsp_par_catalog, dsp_par_cat_file)

psp_par_cat_file = os.path.join(pars_path(setup), "psp", "validity.jsonl")
if os.path.isfile(psp_par_cat_file):
    os.remove(psp_par_cat_file)
pathlib.Path(os.path.dirname(psp_par_cat_file)).mkdir(parents=True, exist_ok=True)
ds.pars_key_resolve.write_to_jsonl(psp_par_catalog, psp_par_cat_file)


localrules:
    gen_filelist,
    autogen_output,


onstart:
    print("INFO: starting workflow")

    # Make sure some packages are initialized before we begin to avoid race conditions
    for pkg in ["dspeed", "lgdo", "matplotlib"]:
        shell('{swenv} python3 -B -c "import ' + pkg + '"')

        # Log parameter catalogs in validity.jsonl files
    hit_par_cat_file = os.path.join(pars_path(setup), "hit", "validity.jsonl")
    if os.path.isfile(hit_par_cat_file):
        os.remove(os.path.join(pars_path(setup), "hit", "validity.jsonl"))
    pathlib.Path(os.path.dirname(hit_par_cat_file)).mkdir(parents=True, exist_ok=True)
    ds.pars_key_resolve.write_to_jsonl(hit_par_catalog, hit_par_cat_file)

    pht_par_cat_file = os.path.join(pars_path(setup), "pht", "validity.jsonl")
    if os.path.isfile(pht_par_cat_file):
        os.remove(os.path.join(pars_path(setup), "pht", "validity.jsonl"))
    pathlib.Path(os.path.dirname(pht_par_cat_file)).mkdir(parents=True, exist_ok=True)
    ds.pars_key_resolve.write_to_jsonl(pht_par_catalog, pht_par_cat_file)

    dsp_par_cat_file = os.path.join(pars_path(setup), "dsp", "validity.jsonl")
    if os.path.isfile(dsp_par_cat_file):
        os.remove(dsp_par_cat_file)
    pathlib.Path(os.path.dirname(dsp_par_cat_file)).mkdir(parents=True, exist_ok=True)
    ds.pars_key_resolve.write_to_jsonl(dsp_par_catalog, dsp_par_cat_file)

    psp_par_cat_file = os.path.join(pars_path(setup), "psp", "validity.jsonl")
    if os.path.isfile(psp_par_cat_file):
        os.remove(psp_par_cat_file)
    pathlib.Path(os.path.dirname(psp_par_cat_file)).mkdir(parents=True, exist_ok=True)
    ds.pars_key_resolve.write_to_jsonl(psp_par_catalog, psp_par_cat_file)


onsuccess:
    from snakemake.report import auto_report

    rep_dir = f"{log_path(setup)}/report-{datetime.strftime(datetime.utcnow(), '%Y%m%dT%H%M%SZ')}"
    pathlib.Path(rep_dir).mkdir(parents=True, exist_ok=True)
    # auto_report(workflow.persistence.dag, f"{rep_dir}/report.html")

    with open(os.path.join(rep_dir, "dag.txt"), "w") as f:
        f.writelines(str(workflow.persistence.dag))
        # shell(f"cat {rep_dir}/dag.txt | dot -Tpdf > {rep_dir}/dag.pdf")

    with open(f"{rep_dir}/rg.txt", "w") as f:
        f.writelines(str(workflow.persistence.dag.rule_dot()))
        # shell(f"cat {rep_dir}/rg.txt | dot -Tpdf > {rep_dir}/rg.pdf")

        # remove .gen files
    files = glob.glob("*.gen")
    for file in files:
        if os.path.isfile(file):
            os.remove(file)

            # remove filelists
    files = glob.glob(os.path.join(filelist_path(setup), "*"))
    for file in files:
        if os.path.isfile(file):
            os.remove(file)
    if os.path.exists(filelist_path(setup)):
        os.rmdir(filelist_path(setup))

        # remove logs
    files = glob.glob(os.path.join(tmp_log_path(setup), "*", "*.log"))
    for file in files:
        if os.path.isfile(file):
            os.remove(file)
    dirs = glob.glob(os.path.join(tmp_log_path(setup), "*"))
    for d in dirs:
        if os.path.isdir(d):
            os.rmdir(d)
    if os.path.exists(tmp_log_path(setup)):
        os.rmdir(tmp_log_path(setup))


rule gen_filelist:
    """Generate file list.

    It is a checkpoint so when it is run it will update the dag passed on the
    files it finds as an output. It does this by taking in the search pattern,
    using this to find all the files that match this pattern, deriving the keys
    from the files found and generating the list of new files needed.
    """
    input:
        lambda wildcards: get_filelist(
            wildcards,
            setup,
            get_pattern_tier_raw(setup),
            ignore_keys_file=os.path.join(configs, "ignore_keys.keylist"),
            analysis_runs_file=os.path.join(configs, "analysis_runs.json"),
        ),
    output:
        os.path.join(filelist_path(setup), "{label}-{tier}.filelist"),
    run:
        if len(input) == 0:
            print(
                "WARNING: No files found for the given pattern\nmake sure pattern follows the format: all-{experiment}-{period}-{run}-{datatype}-{timestamp}-{tier}.gen"
            )
        with open(output[0], "w") as f:
            for fn in input:
                f.write(f"{fn}\n")
