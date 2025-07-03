"""
Snakemake rules for calibrating daq energy for blinding. Two steps:
- determination of calibration curve for daq energy
- combining all channels into single par file
"""

from legenddataflow.methods.patterns import (
    get_pattern_pars,
    get_pattern_plts,
    get_pattern_pars_tmp_channel,
    get_pattern_plts_tmp_channel,
    get_pattern_log_channel,
)
from pathlib import Path
from legenddataflowscripts.workflow import execenv_pyexe


rule build_blinding_calibration:
    """
    Runs a check on the daqenergy of the calibration run that the blinding curve given still applies,
    if so creates a file whose existence will be checked by the raw blinding before proceeding with blinding the phy data
    """
    input:
        files=Path(filelist_path(config))
        / "all-{experiment}-{period}-{run}-cal-raw.filelist",
    params:
        timestamp="{timestamp}",
        datatype="cal",
        channel="{channel}",
        configs=ro(config_path(config)),
        meta=ro(metadata_path(config)),
        raw_table_name=lambda wildcards: get_table_name(
            channelmap_textdb,
            config,
            "cal",
            wildcards.timestamp,
            wildcards.channel,
            "raw",
        ),
    output:
        par_file=temp(get_pattern_pars_tmp_channel(config, "raw_blindcal")),
        plot_file=temp(get_pattern_plts_tmp_channel(config, "raw_blindcal")),
    log:
        get_pattern_log_channel(config, "pars_hit_blind_cal", time),
    group:
        "par-raw-blinding"
    resources:
        runtime=300,
    shell:
        execenv_pyexe(config, "par-geds-raw-blindcal") + "--log {log} "
        "--datatype {params.datatype} "
        "--timestamp {params.timestamp} "
        "--channel {params.channel} "
        "--configs {params.configs} "
        "--meta {params.meta} "
        "--raw-table-name {params.raw_table_name} "
        "--plot-file {output.plot_file} "
        "--blind-curve {output.par_file} "
        "--files {input.files} "


rule build_plts_blinding:
    input:
        lambda wildcards: get_plt_chanlist(
            config,
            f"all-{wildcards.experiment}-{wildcards.period}-{wildcards.run}-cal-{wildcards.timestamp}-channels",
            "raw",
            det_status_textdb,
            channelmap_textdb,
            system="geds",
            name="blindcal",
        ),
    output:
        get_pattern_plts(config, "raw", name="blindcal"),
    group:
        "merge-blindcal"
    shell:
        execenv_pyexe(config, "merge-channels") + "--input {input} "
        "--output {output} "


rule build_pars_blinding:
    input:
        infiles=lambda wildcards: get_par_chanlist(
            config,
            f"all-{wildcards.experiment}-{wildcards.period}-{wildcards.run}-cal-{wildcards.timestamp}-channels",
            "raw",
            det_status_textdb,
            channelmap_textdb,
            extension="yaml",
            system="geds",
            name="blindcal",
        ),
        plts=get_pattern_plts(config, "raw", name="blindcal"),
    output:
        get_pattern_pars(config, "raw", name="blindcal", check_in_cycle=check_in_cycle),
    group:
        "merge-blindcal"
    shell:
        execenv_pyexe(config, "merge-channels") + "--input {input.infiles} "
        "--output {output} "
