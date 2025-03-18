"""DSP parameter generation for SiPM data"""

from pathlib import Path

from legenddataflow import patterns as patt
from legenddataflow import utils, execenv_pyexe


rule build_pars_dsp_tau_spms:
    input:
        raw_file=get_pattern_tier(config, "raw", check_in_cycle=False),
        pardb=lambda wildcards: get_input_par_file(config, wildcards, "dsp", "par_dsp"),
    params:
        timestamp="{timestamp}",
        datatype="{datatype}",
        channels=lambda wildcards: get_chanlist(
            config,
            wildcards.timestamp,
            wildcards.datatype,
            workflow,
            det_status_textdb,
            channelmap_textdb,
            system="spms",
        ),
        raw_table_names=lambda wildcards: [
            get_table_name(
                channelmap_textdb,
                wildcards.datatype,
                wildcards.timestamp,
                channel,
                "raw",
            )
            for channel in get_chanlist(
                config,
                wildcards.timestamp,
                wildcards.datatype,
                workflow,
                det_status_textdb,
                channelmap_textdb,
                system="spms",
            )
        ],
    wildcard_constraints:
        datatype=r"\b(?!cal\b|xtc\b)\w+\b",
    output:
        patt.get_pattern_pars(config, "dsp", name="spms", datatype="{datatype}"),
    log:
        patt.get_pattern_log(config, "pars_spms", time),
    group:
        "par-dsp"
    shell:
        execenv_pyexe(config, "par-spms-dsp-trg-thr-multi") + "--config-path {configs} "
        "--raw-file {input.raw_file} "
        "--dsp-db {input.pardb} "
        "--datatype {params.datatype} "
        "--timestamp {params.timestamp} "
        "--sipm-names {params.channels} "
        "--raw-table-names {params.raw_table_names} "
        "--output-file {output} "
        "--logfile {log} "
