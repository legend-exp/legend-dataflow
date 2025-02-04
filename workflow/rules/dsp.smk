"""
Snakemake rules for processing dsp tier.
- combining of all channels into single pars files with associated plot and results files
- running dsp over all channels using par file
"""

from legenddataflow.pars_loading import ParsCatalog
from legenddataflow.create_pars_keylist import ParsKeyResolve
from pathlib import Path
from legenddataflow.create_pars_keylist import ParsKeyResolve
from legenddataflow.patterns import (
    get_pattern_plts,
    get_pattern_tier,
    get_pattern_pars_tmp,
    get_pattern_log,
    get_pattern_pars,
)
from legenddataflow.execenv import execenv_smk_py_script

dsp_par_catalog = ParsKeyResolve.get_par_catalog(
    ["-*-*-*-cal"],
    get_pattern_tier(config, "raw", check_in_cycle=False),
    {"cal": ["par_dsp"], "lar": ["par_dsp"]},
)

dsp_par_cat_file = Path(pars_path(config)) / "dsp" / "validity.yaml"
if dsp_par_cat_file.is_file():
    dsp_par_cat_file.unlink()
Path(dsp_par_cat_file).parent.mkdir(parents=True, exist_ok=True)
ParsKeyResolve.write_to_yaml(dsp_par_catalog, dsp_par_cat_file)


rule build_plts_dsp:
    input:
        lambda wildcards: get_plt_chanlist(
            config,
            f"all-{wildcards.experiment}-{wildcards.period}-{wildcards.run}-cal-{wildcards.timestamp}-channels",
            "dsp",
            basedir,
            det_status,
            chan_maps,
        ),
    params:
        timestamp="{timestamp}",
        datatype="cal",
    output:
        get_pattern_plts(config, "dsp"),
    group:
        "merge-dsp"
    shell:
        f'{execenv_smk_py_script(config, "merge_channels")}'
        "--input {input} "
        "--output {output} "
        "--channelmap {meta} "


rule build_pars_dsp_objects:
    input:
        lambda wildcards: get_par_chanlist(
            config,
            f"all-{wildcards.experiment}-{wildcards.period}-{wildcards.run}-cal-{wildcards.timestamp}-channels",
            "dsp",
            basedir,
            det_status,
            chan_maps,
            name="objects",
            extension="pkl",
        ),
    params:
        timestamp="{timestamp}",
        datatype="cal",
    output:
        get_pattern_pars(
            config,
            "dsp",
            name="objects",
            extension="dir",
            check_in_cycle=check_in_cycle,
        ),
    group:
        "merge-dsp"
    shell:
        f'{execenv_smk_py_script(config, "merge_channels")}'
        "--input {input} "
        "--output {output} "
        "--timestamp {params.timestamp} "
        "--channelmap {meta} "


rule build_pars_dsp_db:
    input:
        lambda wildcards: get_par_chanlist(
            config,
            f"all-{wildcards.experiment}-{wildcards.period}-{wildcards.run}-cal-{wildcards.timestamp}-channels",
            "dsp",
            basedir,
            det_status,
            chan_maps,
        ),
    params:
        timestamp="{timestamp}",
        datatype="cal",
    output:
        temp(
            get_pattern_pars_tmp(
                config,
                "dsp",
                datatype="cal",
            )
        ),
    group:
        "merge-dsp"
    shell:
        f'{execenv_smk_py_script(config, "merge_channels")}'
        "--input {input} "
        "--output {output} "
        "--timestamp {params.timestamp} "
        "--channelmap {meta} "


rule build_pars_dsp:
    input:
        in_files=lambda wildcards: get_par_chanlist(
            config,
            f"all-{wildcards.experiment}-{wildcards.period}-{wildcards.run}-cal-{wildcards.timestamp}-channels",
            "dsp",
            basedir,
            det_status,
            chan_maps,
            name="dplms",
            extension="lh5",
        ),
        in_db=get_pattern_pars_tmp(
            config,
            "dsp",
            datatype="cal",
        ),
        plts=get_pattern_plts(config, "dsp"),
        objects=get_pattern_pars(
            config,
            "dsp",
            name="objects",
            extension="dir",
            check_in_cycle=check_in_cycle,
        ),
    params:
        timestamp="{timestamp}",
        datatype="cal",
    output:
        out_file=get_pattern_pars(
            config,
            "dsp",
            extension="lh5",
            check_in_cycle=check_in_cycle,
        ),
        out_db=get_pattern_pars(config, "dsp", check_in_cycle=check_in_cycle),
    group:
        "merge-dsp"
    shell:
        f'{execenv_smk_py_script(config, "merge_channels")}'
        "--output {output.out_file} "
        "--in_db {input.in_db} "
        "--out_db {output.out_db} "
        "--input {input.in_files} "
        "--timestamp {params.timestamp} "
        "--channelmap {meta} "


rule build_dsp:
    input:
        raw_file=get_pattern_tier(config, "raw", check_in_cycle=False),
        pars_file=ancient(
            lambda wildcards: ParsCatalog.get_par_file(
                config, wildcards.timestamp, "dsp"
            )
        ),
    params:
        timestamp="{timestamp}",
        datatype="{datatype}",
        ro_input=lambda _, input: {k: ro(v) for k, v in input.items()},
    output:
        tier_file=get_pattern_tier(config, "dsp", check_in_cycle=check_in_cycle),
        db_file=get_pattern_pars_tmp(config, "dsp_db"),
    log:
        get_pattern_log(config, "tier_dsp", time),
    group:
        "tier-dsp"
    resources:
        runtime=300,
        mem_swap=lambda wildcards: 35 if wildcards.datatype == "cal" else 25,
    shell:
        f'{execenv_smk_py_script(config, "build_tier_dsp")}'
        "--log {log} "
        "--tier dsp "
        f"--configs {ro(configs)} "
        "--metadata {meta} "
        "--datatype {params.datatype} "
        "--timestamp {params.timestamp} "
        "--input {params.ro_input[raw_file]} "
        "--output {output.tier_file} "
        "--db_file {output.db_file} "
        "--pars_file {params.ro_input[pars_file]} "
