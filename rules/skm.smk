"""
Snakemake rules for processing skm tier.
"""

from scripts.util.patterns import (
    get_pattern_tier,
    get_pattern_log,
    get_pattern_pars,
)


rule build_skm:
    input:
        dsp_files=os.path.join(
            filelist_path(setup), "all-{experiment}-{period}-{run}-phy-dsp.filelist"
        ),
        hit_files=os.path.join(
            filelist_path(setup), "all-{experiment}-{period}-{run}-phy-pht.filelist"
        ),
        tcm_files=os.path.join(
            filelist_path(setup), "all-{experiment}-{period}-{run}-phy-tcm.filelist"
        ),
        evt_files=lambda wildcards: read_filelist_phy(wildcards, "pet"),
    output:
        skm_file=get_pattern_tier(setup, "skm", check_in_cycle=check_in_cycle),
    params:
        timestamp="{timestamp}",
        datatype="{datatype}",
    log:
        get_pattern_log(setup, "tier_skm"),
    group:
        "tier-skm"
    resources:
        runtime=300,
    shell:
        "{swenv} python3 -B "
        f"{basedir}/../scripts/build_skm.py "
        "--configs {configs} "
        "--metadata {meta} "
        "--log {log} "
        "--datatype {params.datatype} "
        "--timestamp {params.timestamp} "
        "--hit_files {input.hit_files} "
        "--tcm_files {input.tcm_files} "
        "--dsp_files {input.dsp_files} "
        "--evt_files {input.evt_files} "
        "--output {output.skm_file} "
