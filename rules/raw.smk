from scripts.util.patterns import (
    get_pattern_tier_daq,
    get_pattern_tier_raw,
    get_pattern_tier,
    get_pattern_log,
    get_pattern_tier_raw_blind,
)


rule build_raw:
    """
    This rule runs build raw, it takes in a daq file and outputs a raw file
    """
    input:
        get_pattern_tier_daq(setup),
    params:
        timestamp="{timestamp}",
        datatype="{datatype}",
        ro_input=lambda _, input: ro(input),
    output:
        get_pattern_tier(setup, "raw", check_in_cycle=check_in_cycle),
    log:
        get_pattern_log(setup, "tier_raw"),
    group:
        "tier-raw"
    resources:
        mem_swap=110,
        runtime=300,
    shell:
        "{swenv} python3 -B "
        "{basedir}/../scripts/build_raw.py "
        "--log {log} "
        f"--configs {ro(configs)} "
        f"--chan_maps {ro(chan_maps)} "
        "--datatype {params.datatype} "
        "--timestamp {params.timestamp} "
        "{params.ro_input} {output}"


rule build_raw_blind:
    """
    This rule runs the data blinding, it takes in the raw file, calibration curve stored in the overrides
    and runs only if the blinding check file is on disk. Output is just the blinded raw file.
    """
    input:
        tier_file=get_pattern_tier(setup, "raw", check_in_cycle=False).replace(
            "{datatype}", "phy"
        ),
        blind_file=get_blinding_curve_file,
    params:
        timestamp="{timestamp}",
        datatype="phy",
        ro_input=lambda _, input: {k: ro(v) for k, v in input.items()},
    output:
        get_pattern_tier_raw_blind(setup),
    log:
        get_pattern_log(setup, "tier_raw_blind").replace("{datatype}", "phy"),
    group:
        "tier-raw"
    resources:
        mem_swap=110,
        runtime=300,
    shell:
        "{swenv} python3 -B "
        "{basedir}/../scripts/build_raw_blind.py "
        "--log {log} "
        f"--configs {ro(configs)} "
        f"--chan_maps {ro(chan_maps)} "
        f"--metadata {ro(meta)} "
        "--datatype {params.datatype} "
        "--timestamp {params.timestamp} "
        "--blind_curve {params.ro_input[blind_file]} "
        "--input {params.ro_input[tier_file]} "
        "--output {output}"
