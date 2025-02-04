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
from legenddataflow.execenv import execenv_smk_py_script

psp_par_catalog = ParsKeyResolve.get_par_catalog(
    ["-*-*-*-cal"],
    get_pattern_tier(config, "raw", check_in_cycle=False),
    {"cal": ["par_psp"], "lar": ["par_psp"]},
)

psp_par_cat_file = Path(pars_path(config)) / "psp" / "validity.yaml"
if psp_par_cat_file.is_file():
    psp_par_cat_file.unlink()
Path(psp_par_cat_file).parent.mkdir(parents=True, exist_ok=True)
ParsKeyResolve.write_to_yaml(psp_par_catalog, psp_par_cat_file)


include: "channel_merge.smk"


build_merge_rules("psp", lh5_merge=True)


rule build_psp:
    input:
        raw_file=get_pattern_tier(config, "raw", check_in_cycle=False),
        pars_file=ancient(
            lambda wildcards: ParsCatalog.get_par_file(
                config, wildcards.timestamp, "psp"
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
        f'{execenv_smk_py_script(config, "build_tier_dsp")}'
        "--log {log} "
        "--tier psp "
        f"--configs {ro(configs)} "
        "--metadata {meta} "
        "--datatype {params.datatype} "
        "--timestamp {params.timestamp} "
        "--input {params.ro_input[raw_file]} "
        "--output {output.tier_file} "
        "--db_file {output.db_file} "
        "--pars_file {params.ro_input[pars_file]} "
