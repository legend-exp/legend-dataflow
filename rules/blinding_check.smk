"""
Snakemake rules for checking blinding. Two steps:
- check for each channel that 583 and 2614 peaks within 5keV if so produce check file
- combining all channel check files into single check file
"""


def get_blinding_curve(wildcards):
    # func to load in blinding curves
    par_files = ds.pars_catalog.get_calib_files(
        Path(par_overwrite_path(setup)) / "raw" / "validity.jsonl", wildcards.timestamp
    )
    if isinstance(par_files, str):
        return str(Path(par_overwrite_path(setup)) / "raw" / par_files)
    else:
        return [
            str(Path(par_overwrite_path(setup)) / "raw" / par_file)
            for par_file in par_files
        ]


rule build_blinding_check:
    """
    Runs a check on the daqenergy of the calibration run that the blinding curve given still applies,
    if so creates a file whose existence will be checked by the raw blinding before proceeding with blinding the phy data
    """
    input:
        files=lambda wildcards: read_filelist_cal(wildcards, "dsp"),
        par_file=get_blinding_curve,
    params:
        timestamp="{timestamp}",
        datatype="cal",
        channel="{channel}",
    output:
        check_file=temp(get_pattern_pars_tmp_channel(setup, "raw")),
        plot_file=temp(get_pattern_plts_tmp_channel(setup, "raw")),
    log:
        get_pattern_log_channel(setup, "pars_hit_blind_check"),
    group:
        "par-hit"
    resources:
        runtime=300,
    shell:
        "{swenv} python3 -B "
        f"{workflow.source_path('../scripts/check_blinding.py')} "
        "--log {log} "
        "--datatype {params.datatype} "
        "--timestamp {params.timestamp} "
        "--channel {params.channel} "
        "--configs {configs} "
        "--output {output.check_file} "
        "--blind_curve {input.par_file} "
        "--plot_file {output.plot_file} "
        "--files {input.files} "


checkpoint build_pars_raw:
    input:
        lambda wildcards: read_filelist_pars_cal_channel(wildcards, "raw"),
        lambda wildcards: read_filelist_plts_cal_channel(wildcards, "raw"),
    output:
        get_pattern_par_raw(setup),
        get_pattern_plts(setup, "raw"),
    group:
        "merge-blinding"
    shell:
        "{swenv} python3 -B "
        f"{workflow.source_path('../scripts/merge_channels.py')} "
        "--input {input} "
        "--output {output} "
