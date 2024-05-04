"""
Snakemake file containing the rules for generating the tcm
"""

from scripts.util.patterns import (
    get_pattern_tier_raw,
    get_pattern_tier,
    get_pattern_log,
    get_pattern_pars_tmp_channel,
    get_pattern_log_channel,
)


# This rule builds the tcm files each raw file
rule build_tier_tcm:
    input:
        get_pattern_tier_raw(setup),
    params:
        timestamp="{timestamp}",
        datatype="{datatype}",
    output:
        get_pattern_tier(setup, "tcm", check_in_cycle=check_in_cycle),
    log:
        get_pattern_log(setup, "tier_tcm"),
    group:
        "tier-tcm"
    resources:
        runtime=300,
    shell:
        "{swenv} python3 -B "
        "{basedir}/../scripts/build_tcm.py "
        "--log {log} "
        "--configs {configs} "
        "--datatype {params.datatype} "
        "--timestamp {params.timestamp} "
        "{input} "
        "{output}"


# This rule builds the tcm files each raw file
rule build_pulser_ids:
    input:
        tcm_files=os.path.join(
            filelist_path(setup), "all-{experiment}-{period}-{run}-cal-tcm.filelist"
        ),
    params:
        timestamp="{timestamp}",
        datatype="cal",
        channel="{channel}",
    output:
        pulser=temp(get_pattern_pars_tmp_channel(setup, "tcm", "pulser_ids")),
    log:
        get_pattern_log_channel(setup, "tcm_pulsers"),
    group:
        "tier-tcm"
    resources:
        runtime=300,
    shell:
        "{swenv} python3 -B "
        "{basedir}/../scripts/pars_tcm_pulser.py "
        "--log {log} "
        "--configs {configs} "
        "--datatype {params.datatype} "
        "--timestamp {params.timestamp} "
        "--channel {params.channel} "
        "--tcm_files {input.tcm_files} "
        "--pulser_file {output.pulser} "
