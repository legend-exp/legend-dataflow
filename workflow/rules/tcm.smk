"""
Snakemake file containing the rules for generating the tcm
"""

from legenddataflow.patterns import (
    get_pattern_tier,
    get_pattern_log,
    get_pattern_pars_tmp_channel,
    get_pattern_log_channel,
)
from legenddataflow.execenv import execenv_pyexe


# This rule builds the tcm files each raw file
rule build_tier_tcm:
    input:
        get_pattern_tier(config, "raw", check_in_cycle=False),
    params:
        timestamp="{timestamp}",
        datatype="{datatype}",
        input=lambda _, input: ro(input),
    output:
        get_pattern_tier(config, "tcm", check_in_cycle=check_in_cycle),
    log:
        get_pattern_log(config, "tier_tcm", time),
    group:
        "tier-tcm"
    resources:
        runtime=300,
        mem_swap=20,
    shell:
        execenv_pyexe(config, "build-tier-tcm") + "--log {log} "
        f"--configs {ro(configs)} "
        "--datatype {params.datatype} "
        "--timestamp {params.timestamp} "
        "-- {params.input} {output}"


# This rule builds the tcm files each raw file
rule build_pulser_ids:
    input:
        os.path.join(
            filelist_path(config), "all-{experiment}-{period}-{run}-cal-tcm.filelist"
        ),
    params:
        input=lambda _, input: ro(input),
        timestamp="{timestamp}",
        datatype="cal",
        channel="{channel}",
        rawid=lambda wildcards: channelmap_textdb.valid_for(
            wildcards.timestamp, system="cal"
        )[wildcards.channel].daq.rawid,
    output:
        pulser=temp(get_pattern_pars_tmp_channel(config, "tcm", "pulser_ids")),
    log:
        get_pattern_log_channel(config, "tcm_pulsers", time),
    group:
        "tier-tcm"
    resources:
        runtime=300,
    shell:
        execenv_pyexe(config, "par-geds-tcm-pulser") + "--log {log} "
        f"--configs {ro(configs)} "
        "--datatype {params.datatype} "
        "--timestamp {params.timestamp} "
        "--channel {params.channel} "
        "--rawid {params.rawid} "
        "--tcm-files {params.input} "
        "--pulser-file {output.pulser} "
        "--metadata {meta} "
