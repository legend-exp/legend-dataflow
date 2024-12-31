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

from pathlib import Path
import os
import json
import sys
import glob
from datetime import datetime
from collections import OrderedDict
import logging

import scripts.util as ds
from scripts.util.pars_loading import ParsCatalog
from scripts.util.patterns import get_pattern_tier
from scripts.util.utils import (
    subst_vars_in_snakemake_config,
    runcmd,
    config_path,
    chan_map_path,
    filelist_path,
    metadata_path,
    tmp_log_path,
    pars_path,
    det_status_path,
)

# Set with `snakemake --configfile=/path/to/your/config.json`
# configfile: "have/to/specify/path/to/your/config.json"

subst_vars_in_snakemake_config(workflow, config)

check_in_cycle = True
setup = config["setups"]["l200"]
configs = config_path(setup)
chan_maps = chan_map_path(setup)
meta = metadata_path(setup)
det_status = det_status_path(setup)
swenv = runcmd(setup)
part = ds.CalGrouping(setup, Path(det_status) / "cal_groupings.yaml")
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
include: "rules/ann.smk"
include: "rules/evt.smk"
include: "rules/skm.smk"
include: "rules/blinding_calibration.smk"
include: "rules/qc_phy.smk"


localrules:
    gen_filelist,
    autogen_output,


onstart:
    print("INFO: starting workflow")

    # Make sure some packages are initialized before we begin to avoid race conditions
    for pkg in ["dspeed", "lgdo", "matplotlib"]:
        shell('{swenv} python3 -B -c "import ' + pkg + '"')

        # Log parameter catalogs in validity.jsonl files
    hit_par_cat_file = Path(pars_path(setup)) / "hit" / "validity.yaml"
    if hit_par_cat_file.is_file():
        hit_par_cat_file.unlink()
    Path(hit_par_cat_file).parent.mkdir(parents=True, exist_ok=True)
    ds.ParsKeyResolve.write_to_yaml(hit_par_catalog, hit_par_cat_file)

    pht_par_cat_file = Path(pars_path(setup)) / "pht" / "validity.yaml"
    if pht_par_cat_file.is_file():
        pht_par_cat_file.unlink()
    Path(pht_par_cat_file).parent.mkdir(parents=True, exist_ok=True)
    ds.ParsKeyResolve.write_to_yaml(pht_par_catalog, pht_par_cat_file)

    dsp_par_cat_file = Path(pars_path(setup)) / "dsp" / "validity.yaml"
    if dsp_par_cat_file.is_file():
        dsp_par_cat_file.unlink()
    Path(dsp_par_cat_file).parent.mkdir(parents=True, exist_ok=True)
    ds.ParsKeyResolve.write_to_yaml(dsp_par_catalog, dsp_par_cat_file)

    psp_par_cat_file = Path(pars_path(setup)) / "psp" / "validity.yaml"
    if psp_par_cat_file.is_file():
        psp_par_cat_file.unlink()
    Path(psp_par_cat_file).parent.mkdir(parents=True, exist_ok=True)
    ds.ParsKeyResolve.write_to_yaml(psp_par_catalog, psp_par_cat_file)


onsuccess:
    from snakemake.report import auto_report

    rep_dir = f"{log_path(setup)}/report-{datetime.strftime(datetime.utcnow(), '%Y%m%dT%H%M%SZ')}"
    Path(rep_dir).mkdir(parents=True, exist_ok=True)
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

            #  remove filelists
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
            get_pattern_tier(setup, "raw", check_in_cycle=False),
            ignore_keys_file=Path(det_status) / "ignored_daq_cycles.yaml",
            analysis_runs_file=Path(det_status) / "runlists.yaml",
        ),
    output:
        temp(Path(filelist_path(setup)) / "{label}-{tier}.filelist"),
    run:
        if len(input) == 0:
            print(
                f"WARNING: No files found for the given pattern:{wildcards.label}",
                "\nmake sure pattern follows the format: all-{experiment}-{period}-{run}-{datatype}-{timestamp}-{tier}.gen",
            )
        with open(output[0], "w") as f:
            for fn in input:
                f.write(f"{fn}\n")
