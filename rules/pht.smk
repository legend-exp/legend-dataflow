"""
Snakemake rules for processing pht (partition hit) tier data. This is done in 4 steps:
- extraction of calibration curves(s) for each run for each channel from cal data
- extraction of psd calibration parameters and partiton level energy fitting for each channel over whole partition from cal data 
- combining of all channels into single pars files with associated plot and results files
- running build hit over all channels using par file
"""

# This rule builds the energy calibration using the calibration dsp files
rule build_per_energy_calibration:
    input:
        files=lambda wildcards: read_filelist_cal(wildcards, "dsp"),
        ctc_dict=ancient(lambda wildcards : ds.pars_catalog.get_par_file(setup, wildcards.timestamp, "dsp")),
    params:
        timestamp="{timestamp}",
        datatype="cal",
        channel="{channel}",
    output:
        ecal_file=temp(get_pattern_pars_tmp_channel(setup, "pht", "energy_cal")),
        results_file=temp(
            get_pattern_pars_tmp_channel(
                setup, "pht", "energy_cal_results", extension="pkl"
            )
        ),
        plot_file=temp(get_pattern_plts_tmp_channel(setup, "pht", "energy_cal")),
    log:
        get_pattern_log_channel(setup, "pars_pht_energy_cal"),
    group:
        "par-pht"
    resources:
        runtime=300,
    shell:
        "{swenv} python3 -B "
        f"{workflow.source_path('../scripts/pars_hit_ecal.py')} "
        "--log {log} "
        "--datatype {params.datatype} "
        "--timestamp {params.timestamp} "
        "--channel {params.channel} "
        "--configs {configs} "
        "--plot_path {output.plot_file} "
        "--results_path {output.results_file} "
        "--save_path {output.ecal_file} "
        "--ctc_dict {input.ctc_dict} "
        "--files {input.files}"


checkpoint build_pars_pht:
    input:
        lambda wildcards: read_filelist_pars_cal_channel(wildcards, 'pht'),
        lambda wildcards: read_filelist_plts_cal_channel(wildcards, 'pht'),
        lambda wildcards: read_filelist_pars_cal_channel(wildcards, 'pht_results'),
    output:
        get_pattern_par_pht(setup),
        get_pattern_par_pht(setup, name="results", extension="dir"),
        get_pattern_plts(setup, "pht"),
    group:
        "merge-hit"
    shell:
        "{swenv} python3 -B "
        f"{workflow.source_path('../scripts/merge_channels.py')} "
        "--input {input} "
        "--output {output} "

rule build_pht:
    input:
        dsp_file=get_pattern_tier_dsp(setup),
        #hit_file = get_pattern_tier_hit(setup),
        pars_file=lambda wildcards : ds.pars_catalog.get_par_file(setup, wildcards.timestamp, "pht"),
    output:
        tier_file=get_pattern_tier_pht(setup),
        db_file=get_pattern_pars_tmp(setup, "pht_db"),
    params:
        timestamp="{timestamp}",
        datatype="{datatype}",
        tier="pht",
    log:
        get_pattern_log(setup, "tier_pht"),
    group:
        "tier-pht"
    resources:
        runtime=300,
    shell:
        "{swenv} python3 -B "
        f"{workflow.source_path('../scripts/build_hit.py')} "
        "--configs {configs} "
        "--log {log} "
        "--tier {params.tier} "
        "--datatype {params.datatype} "
        "--timestamp {params.timestamp} "
        "--pars_file {input.pars_file} "
        "--output {output.tier_file} "
        "--input {input.dsp_file} "
        "--db_file {output.db_file}"


# def fix_name(new_name):
#     """ sets the name of the most recently created rule to be new_name
#     """
#     list(workflow.rules)[-1].name = new_name
#     temp_rules = list(rules.__dict__.items())
#     temp_rules[-1] = (new_name, temp_rules[-1][1])
#     rules.__dict__ = dict(temp_rules)

part_pht_rules = {}
for key, dataset in part.datasets.items():
    for partition in dataset.keys():

        rule:
            input:
                files=part.get_filelists(partition, key, "dsp"),
                ecal_file=part.get_par_files(
                    f"{par_pht_path(setup)}/validity.jsonl",
                    partition,
                    key,
                    tier="pht",
                    name="energy_cal",
                ),
                eres_file=part.get_par_files(
                    f"{par_pht_path(setup)}/validity.jsonl",
                    partition,
                    key,
                    tier="pht",
                    name="energy_cal_results",
                    extension="pkl",
                ),
                inplots=part.get_plt_files(
                    f"{par_pht_path(setup)}/validity.jsonl",
                    partition,
                    key,
                    tier="pht",
                    name="energy_cal",
                ),
            params:
                datatype="cal",
                channel="{channel}",
                timestamp=part.get_timestamp(
                    f"{par_pht_path(setup)}/validity.jsonl", partition, key, tier="pht"
                ),
            output:
                hit_pars=[
                    temp(file)
                    for file in part.get_par_files(
                        f"{par_pht_path(setup)}/validity.jsonl",
                        partition,
                        key,
                        tier="pht",
                    )
                ],
                aoe_results=[
                    temp(file)
                    for file in part.get_par_files(
                        f"{par_pht_path(setup)}/validity.jsonl",
                        partition,
                        key,
                        tier="pht",
                        name="results",
                        extension="pkl",
                    )
                ],
                plot_file=[
                    temp(file)
                    for file in part.get_plt_files(
                        f"{par_pht_path(setup)}/validity.jsonl",
                        partition,
                        key,
                        tier="pht",
                    )
                ],
            log:
                part.get_log_file(
                    f"{par_pht_path(setup)}/validity.jsonl",
                    partition,
                    key,
                    "pht",
                    name="par_pht",
                ),
            group:
                "par-pht"
            resources:
                mem_swap=75,
                runtime=300,
            shell:
                "{swenv} python3 -B "
                f"{workflow.source_path('../scripts/pars_pht.py')} "
                "--log {log} "
                "--configs {configs} "
                "--datatype {params.datatype} "
                "--timestamp {params.timestamp} "
                "--inplots {input.inplots} "
                "--channel {params.channel} "
                "--aoe_results {output.aoe_results} "
                "--eres_file {input.eres_file} "
                "--hit_pars {output.hit_pars} "
                "--plot_file {output.plot_file} "
                "--ecal_file {input.ecal_file} "
                "--input_files {input.files}"

        # fix_name(f"{key}-{partition}")

        if key in part_pht_rules:
            part_pht_rules[key].append(list(workflow.rules)[-1])
        else:
            part_pht_rules[key] = [list(workflow.rules)[-1]]


# Merged energy and a/e supercalibrations to reduce number of rules as they have same inputs/outputs
# This rule builds the a/e calibration using the calibration dsp files for the whole partition
rule build_pht_super_calibrations:
    input:
        files=os.path.join(
            filelist_path(setup), "all-{experiment}-{period}-{run}-cal-dsp.filelist"
        ),
        ecal_file=get_pattern_pars_tmp_channel(setup, "pht", "energy_cal"),
        eres_file=get_pattern_pars_tmp_channel(
            setup, "pht", "energy_cal_results", extension="pkl"
        ),
        inplots=get_pattern_plts_tmp_channel(setup, "pht", "energy_cal"),
    params: 
        datatype="cal",
        channel="{channel}",
        timestamp="{timestamp}",
    output:
        hit_pars=temp(get_pattern_pars_tmp_channel(setup, "pht")),
        aoe_results=temp(
            get_pattern_pars_tmp_channel(setup, "pht", "results", extension="pkl")
        ),
        plot_file=temp(get_pattern_plts_tmp_channel(setup, "pht")),
    log:
        get_pattern_log_channel(setup, "pars_pht_aoe_cal"),
    group:
        "par-pht"
    resources:
        mem_swap=60,
        runtime=300,
    shell:
        "{swenv} python3 -B "
        f"{workflow.source_path('../scripts/pars_pht.py')} "
        "--log {log} "
        "--configs {configs} "
        "--datatype {params.datatype} "
        "--timestamp {params.timestamp} "
        "--inplots {input.inplots} "
        "--channel {params.channel} "
        "--aoe_results {output.aoe_results} "
        "--eres_file {input.eres_file} "
        "--hit_pars {output.hit_pars} "
        "--plot_file {output.plot_file} "
        "--ecal_file {input.ecal_file} "
        "--input_files {input.files}"


fallback_pht_rule = list(workflow.rules)[-1]

rule_order_list = []
ordered = OrderedDict(part_pht_rules)
ordered.move_to_end("default")
for key, items in ordered.items():
    rule_order_list += [item.name for item in items]
rule_order_list.append(fallback_pht_rule.name)
workflow._ruleorder.add(*rule_order_list)  # [::-1]

