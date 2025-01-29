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

psp_par_catalog = ParsKeyResolve.get_par_catalog(
    ["-*-*-*-cal"],
    get_pattern_tier(setup, "raw", check_in_cycle=False),
    {"cal": ["par_psp"], "lar": ["par_psp"]},
)

psp_par_cat_file = Path(pars_path(setup)) / "psp" / "validity.yaml"
if psp_par_cat_file.is_file():
    psp_par_cat_file.unlink()
Path(psp_par_cat_file).parent.mkdir(parents=True, exist_ok=True)
ParsKeyResolve.write_to_yaml(psp_par_catalog, psp_par_cat_file)


rule build_pars_psp_objects:
    input:
        lambda wildcards: get_par_chanlist(
            setup,
            f"all-{wildcards.experiment}-{wildcards.period}-{wildcards.run}-cal-{wildcards.timestamp}-channels",
            "psp",
            basedir,
            det_status,
            chan_maps,
            name="objects",
            extension="pkl",
        ),
    output:
        get_pattern_pars(
            setup,
            "psp",
            name="objects",
            extension="dir",
            check_in_cycle=check_in_cycle,
        ),
    group:
        "merge-psp"
    shell:
        "{swenv} python3 -B "
        "{basedir}/../scripts/merge_channels.py "
        "--input {input} "
        "--output {output} "
        "--channelmap {meta} "


rule build_plts_psp:
    input:
        lambda wildcards: get_plt_chanlist(
            setup,
            f"all-{wildcards.experiment}-{wildcards.period}-{wildcards.run}-cal-{wildcards.timestamp}-channels",
            "psp",
            basedir,
            det_status,
            chan_maps,
        ),
    output:
        get_pattern_plts(setup, "psp"),
    group:
        "merge-psp"
    shell:
        "{swenv} python3 -B "
        "{basedir}/../scripts/merge_channels.py "
        "--input {input} "
        "--output {output} "
        "--channelmap {meta} "


rule build_pars_psp_db:
    input:
        lambda wildcards: get_par_chanlist(
            setup,
            f"all-{wildcards.experiment}-{wildcards.period}-{wildcards.run}-cal-{wildcards.timestamp}-channels",
            "psp",
            basedir,
            det_status,
            chan_maps,
        ),
    output:
        temp(
            get_pattern_pars_tmp(
                setup,
                "psp",
                datatype="cal",
            )
        ),
    group:
        "merge-psp"
    shell:
        "{swenv} python3 -B "
        "{basedir}/../scripts/merge_channels.py "
        "--input {input} "
        "--output {output} "
        "--channelmap {meta} "


rule build_pars_psp:
    input:
        in_files=lambda wildcards: get_par_chanlist(
            setup,
            f"all-{wildcards.experiment}-{wildcards.period}-{wildcards.run}-cal-{wildcards.timestamp}-channels",
            "dsp",
            basedir,
            det_status,
            chan_maps,
            name="dplms",
            extension="lh5",
        ),
        in_db=get_pattern_pars_tmp(
            setup,
            "psp",
            datatype="cal",
        ),
        plts=get_pattern_plts(setup, "psp"),
        objects=get_pattern_pars(
            setup,
            "psp",
            name="objects",
            extension="dir",
            check_in_cycle=check_in_cycle,
        ),
    output:
        out_file=get_pattern_pars(
            setup,
            "psp",
            extension="lh5",
            check_in_cycle=check_in_cycle,
        ),
        out_db=get_pattern_pars(setup, "psp", check_in_cycle=check_in_cycle),
    group:
        "merge-psp"
    shell:
        "{swenv} python3 -B "
        "{basedir}/../scripts/merge_channels.py "
        "--output {output.out_file} "
        "--in_db {input.in_db} "
        "--out_db {output.out_db} "
        "--input {input.in_files} "
        "--channelmap {meta} "


rule build_psp:
    input:
        raw_file=get_pattern_tier(setup, "raw", check_in_cycle=False),
        pars_file=ancient(
            lambda wildcards: ParsCatalog.get_par_file(
                setup, wildcards.timestamp, "psp"
            )
        ),
    params:
        timestamp="{timestamp}",
        datatype="{datatype}",
        ro_input=lambda _, input: {k: ro(v) for k, v in input.items()},
    output:
        tier_file=get_pattern_tier(setup, "psp", check_in_cycle=check_in_cycle),
        db_file=get_pattern_pars_tmp(setup, "psp_db"),
    log:
        get_pattern_log(setup, "tier_psp"),
    group:
        "tier-dsp"
    resources:
        runtime=300,
        mem_swap=lambda wildcards: 35 if wildcards.datatype == "cal" else 25,
    shell:
        "{swenv} python3 -B "
        "{basedir}/../scripts/build_dsp.py "
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
