"""
Snakemake file containing the rules for generating the tcm
"""

from legenddataflow.patterns import (
    get_pattern_tier,
    get_pattern_log,
    get_pattern_pars_tmp_channel,
    get_pattern_log_channel,
)


# This rule builds the tcm files each raw file
rule build_tier_tcm:
    input:
        get_pattern_tier(setup, "raw", check_in_cycle=False),
    params:
        timestamp="{timestamp}",
        datatype="{datatype}",
        input=lambda _, input: ro(input),
    output:
        get_pattern_tier(setup, "tcm", check_in_cycle=check_in_cycle),
    log:
        get_pattern_log(setup, "tier_tcm", time),
    group:
        "tier-tcm"
    resources:
        runtime=300,
        mem_swap=20,
    shell:
        "{swenv} python3 -B "
        "{basedir}/../scripts/build_tcm.py "
        "--log {log} "
        f"--configs {ro(configs)} "
        "--datatype {params.datatype} "
        "--timestamp {params.timestamp} "
        "-- {params.input} {output}"


# This rule builds the tcm files each raw file
rule build_pulser_ids:
    input:
        os.path.join(
            filelist_path(setup), "all-{experiment}-{period}-{run}-cal-tcm.filelist"
        ),
    params:
        input=lambda _, input: ro(input),
        timestamp="{timestamp}",
        datatype="cal",
        channel="{channel}",
    output:
        pulser=temp(get_pattern_pars_tmp_channel(setup, "tcm", "pulser_ids")),
    log:
        get_pattern_log_channel(setup, "tcm_pulsers", time),
    group:
        "tier-tcm"
    resources:
        runtime=300,
    shell:
        "{swenv} python3 -B "
        "{basedir}/../scripts/pars_tcm_pulser.py "
        "--log {log} "
        f"--configs {ro(configs)} "
        "--datatype {params.datatype} "
        "--timestamp {params.timestamp} "
        "--channel {params.channel} "
        "--tcm_files {params.input} "
        "--pulser_file {output.pulser} "
        "--metadata {meta} "
