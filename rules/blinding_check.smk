"""
Snakemake rules for checking blinding. Two steps:
- check for each channel that 583 and 2614 peaks within 5keV if so produce check file
- combining all channel check files into single check file
"""

from scripts.util.patterns import (
    get_pattern_pars_tmp_channel,
    get_pattern_plts_tmp_channel,
    get_pattern_log_channel,
    get_pattern_par_raw,
    get_pattern_plts,
    get_pattern_pars,
)


rule build_blinding_check:
    """
    Runs a check on the daqenergy of the calibration run that the blinding curve given still applies,
    if so creates a file whose existence will be checked by the raw blinding before proceeding with blinding the phy data
    """
    input:
        files=os.path.join(
            filelist_path(setup), "all-{experiment}-{period}-{run}-cal-raw.filelist"
        ),
        par_file=get_blinding_curve_file,
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
        "{basedir}/../scripts/check_blinding.py "
        "--log {log} "
        "--datatype {params.datatype} "
        "--timestamp {params.timestamp} "
        "--channel {params.channel} "
        "--configs {configs} "
        "--metadata {meta} "
        "--output {output.check_file} "
        "--blind_curve {input.par_file} "
        "--plot_file {output.plot_file} "
        "--files {input.files} "


rule build_plts_raw:
    input:
        lambda wildcards: get_plt_chanlist(
            setup,
            f"all-{wildcards.experiment}-{wildcards.period}-{wildcards.run}-cal-{wildcards.timestamp}-channels",
            "raw",
            basedir,
            configs,
            chan_maps,
        ),
    output:
        get_pattern_plts(setup, "raw"),
    group:
        "merge-raw"
    shell:
        "{swenv} python3 -B "
        "{basedir}/../scripts/merge_channels.py "
        "--input {input} "
        "--output {output} "


rule build_pars_raw:
    input:
        infiles=lambda wildcards: get_par_chanlist(
            setup,
            f"all-{wildcards.experiment}-{wildcards.period}-{wildcards.run}-cal-{wildcards.timestamp}-channels",
            "raw",
            basedir,
            configs,
            chan_maps,
        ),
        plts=get_pattern_plts(
            setup,
            "raw",
        ),
    output:
        get_pattern_pars(setup, "raw", check_in_cycle=check_in_cycle),
    group:
        "merge-raw"
    shell:
        "{swenv} python3 -B "
        "{basedir}/../scripts/merge_channels.py "
        "--input {input.infiles} "
        "--output {output} "
