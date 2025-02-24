"""DSP parameter generation for SiPM data"""

from pathlib import Path

from legenddataflow import patterns as patt
from legenddataflow import utils, execenv_pyexe


def _get_first_raw_file(wildcards):
    w = wildcards
    flist = (
        Path(utils.filelist_path(config))
        / f"all-{w.experiment}-{w.period}-{w.run}-{w.datatype}-raw.filelist"
    )
    with flist.open() as f:
        return f.read().splitlines()[0]


rule build_pars_dsp_tau_spms:
    input:
        raw_file=lambda wildcards: _get_first_raw_file(wildcards),
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
        "--raw-file {input.raw_file} "
        "--dsp-db {input.pardb} "
        "--datatype {params.datatype} "
        "--timestamp {params.timestamp} "
        "--sipm-name {params.channel} "
        "--raw-table-name {params.raw_table_name} "
        "--output-file {output} "
        "--logfile {log} "
