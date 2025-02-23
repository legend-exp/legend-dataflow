"""DSP parameter generation for SiPM data"""

from pathlib import Path

from legenddataflow import patterns as patt
from legenddataflow import utils, execenv_pyexe


rule build_pars_dsp_tau_spms:
    input:
        filelist=Path(utils.filelist_path(config))
        / "all-{experiment}-{period}-{run}-{datatype}-raw.filelist",
        pardb=lambda wildcards: get_input_par_file(config, wildcards, "dsp", "par_dsp"),
    params:
        timestamp="{timestamp}",
        datatype="{datatype}",
        channel="{channel}",
        raw_table_name=lambda wildcards: get_table_name(
            metadata,
            config,
            wildcards.datatype,
            wildcards.timestamp,
            wildcards.channel,
            "raw",
        ),
    wildcard_constraints:
        datatype=r"\b(?!cal\b|xtc\b)\w+\b",
    output:
        temp(patt.get_pattern_pars_tmp_channel(config, "dsp", datatype="{datatype}")),
    log:
        patt.get_pattern_log_channel(config, "pars_spms", time, datatype="{datatype}"),
    group:
        "par-dsp"
    shell:
        execenv_pyexe(config, "par-spms-dsp-trg-thr") + "--config-path {configs} "
        "--raw-files {input.filelist} "
        "--dsp-db {input.pardb} "
        "--datatype {params.datatype} "
        "--timestamp {params.timestamp} "
        "--sipm-name {params.channel} "
        "--raw-table-name {params.raw_table_name} "
        "--output-file {output} "
        "--logfile {log} "
