"""
Snakemake rules for checking blinding. Two steps:
- check for each channel that 583 and 2614 peaks within 5keV if so produce check file
- combining all channel check files into single check file
"""

from legenddataflow.methods.patterns import (
    get_pattern_pars_tmp_channel,
    get_pattern_plts_tmp_channel,
    get_pattern_log_channel,
    get_pattern_pars,
    get_pattern_plts,
    get_pattern_pars,
)
from legenddataflow.methods.paths import config_path, metadata_path
from legenddataflowscripts.workflow import execenv_pyexe
from legenddataflow.scripts.flow.build_chanlist import (
    get_par_chanlist,
    get_plt_chanlist,
)
from pathlib import Path


rule build_blinding_check:
    """
    Runs a check on the daqenergy of the calibration run that the blinding curve given still applies,
    if so creates a file whose existence will be checked by the raw blinding before proceeding with blinding the phy data
    """
    input:
        files=Path(filelist_path(config))
        / "all-{experiment}-{period}-{run}-cal-raw.filelist",
        par_file=get_blinding_curve_file,
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
        check_file=temp(get_pattern_pars_tmp_channel(config, "raw")),
        plot_file=temp(get_pattern_plts_tmp_channel(config, "raw")),
    log:
        get_pattern_log_channel(config, "pars_hit_blind_check", time),
    group:
        "par-hit"
    resources:
        runtime=300,
    shell:
        execenv_pyexe(config, "par-geds-raw-blindcheck") + "--log {log} "
        "--datatype {params.datatype} "
        "--timestamp {params.timestamp} "
        "--channel {params.channel} "
        "--configs {params.configs} "
        "--metadata {params.meta} "
        "--raw-table-name {params.raw_table_name} "
        "--output {output.check_file} "
        "--blind-curve {input.par_file} "
        "--plot-file {output.plot_file} "
        "--files {input.files} "


rule build_plts_raw:
    input:
        lambda wildcards: get_plt_chanlist(
            config,
            f"all-{wildcards.experiment}-{wildcards.period}-{wildcards.run}-cal-{wildcards.timestamp}-channels",
            "raw",
            det_status_textdb,
            channelmap_textdb,
            system="geds",
        ),
    output:
        get_pattern_plts(config, "raw"),
    group:
        "merge-raw"
    shell:
        execenv_pyexe(config, "merge-channels") + "--input {input} "
        "--output {output} "


rule build_pars_raw:
    input:
        infiles=lambda wildcards: get_par_chanlist(
            config,
            f"all-{wildcards.experiment}-{wildcards.period}-{wildcards.run}-cal-{wildcards.timestamp}-channels",
            "raw",
            det_status_textdb,
            channelmap_textdb,
            extension="yaml",
            system="geds",
        ),
        plts=get_pattern_plts(
            config,
            "raw",
        ),
    output:
        get_pattern_pars(config, "raw", check_in_cycle=check_in_cycle),
    group:
        "merge-raw"
    shell:
        execenv_pyexe(config, "merge-channels") + "--input {input.infiles} "
        "--output {output} "
