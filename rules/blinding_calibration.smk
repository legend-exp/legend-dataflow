"""
Snakemake rules for calibrating daq energy for blinding. Two steps:
- determination of calibration curve for daq energy
- combining all channels into single par file
"""

from scripts.util.patterns import (
    get_pattern_par_raw,
    get_pattern_plts,
    get_pattern_pars_tmp_channel,
    get_pattern_plts_tmp_channel,
    get_pattern_log_channel,
)


rule build_blinding_calibration:
    """
    Runs a check on the daqenergy of the calibration run that the blinding curve given still applies,
    if so creates a file whose existence will be checked by the raw blinding before proceeding with blinding the phy data
    """
    input:
        files=lambda wildcards: read_filelist_cal(wildcards, "raw"),
    params:
        timestamp="{timestamp}",
        datatype="cal",
        channel="{channel}",
        meta=meta,
    output:
        par_file=temp(get_pattern_pars_tmp_channel(setup, "raw_blindcal")),
        plot_file=temp(get_pattern_plts_tmp_channel(setup, "raw_blindcal")),
    log:
        get_pattern_log_channel(setup, "pars_hit_blind_cal"),
    group:
        "par-raw-blinding"
    resources:
        runtime=300,
    shell:
        "{swenv} python3 -B "
        f"{basedir}/../scripts/blinding_calibration.py "
        "--log {log} "
        "--datatype {params.datatype} "
        "--timestamp {params.timestamp} "
        "--channel {params.channel} "
        "--configs {configs} "
        "--meta {params.meta} "
        "--plot_file {output.plot_file} "
        "--blind_curve {output.par_file} "
        "--files {input.files} "


checkpoint build_pars_blinding:
    input:
        lambda wildcards: read_filelist_pars_cal_channel(wildcards, "raw_blindcal"),
        lambda wildcards: read_filelist_plts_cal_channel(wildcards, "raw_blindcal"),
    output:
        get_pattern_par_raw(setup, name="blindcal"),
        get_pattern_plts(setup, "raw", name="blindcal"),
    group:
        "merge-blinding"
    shell:
        "{swenv} python3 -B "
        f"{basedir}/../scripts/merge_channels.py "
        "--input {input} "
        "--output {output} "
