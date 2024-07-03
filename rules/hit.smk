"""
Snakemake rules for processing hit tier. This is done in 4 steps:
- extraction of calibration curves(s) for each channel from cal data
- extraction of psd calibration parameters for each channel from cal data
- combining of all channels into single pars files with associated plot and results files
- running build hit over all channels using par file
"""

from scripts.util.pars_loading import pars_catalog
from scripts.util.patterns import (
    get_pattern_pars_tmp_channel,
    get_pattern_plts_tmp_channel,
    get_pattern_log_channel,
    get_pattern_par_hit,
    get_pattern_plts,
    get_pattern_tier_dsp,
    get_pattern_tier,
    get_pattern_pars_tmp,
    get_pattern_log,
    get_pattern_pars,
)


# This rule builds the qc using the calibration dsp files and fft files
rule build_qc:
    input:
        files=os.path.join(
            filelist_path(setup), "all-{experiment}-{period}-{run}-cal-dsp.filelist"
        ),
        fft_files=os.path.join(
            filelist_path(setup), "all-{experiment}-{period}-{run}-fft-dsp.filelist"
        ),
        pulser=get_pattern_pars_tmp_channel(setup, "tcm", "pulser_ids"),
    params:
        timestamp="{timestamp}",
        datatype="cal",
        channel="{channel}",
    output:
        qc_file=temp(get_pattern_pars_tmp_channel(setup, "hit", "qc")),
        plot_file=temp(get_pattern_plts_tmp_channel(setup, "hit", "qc")),
    log:
        get_pattern_log_channel(setup, "pars_hit_qc"),
    group:
        "par-hit"
    resources:
        runtime=300,
    shell:
        "{swenv} python3 -B "
        "{basedir}/../scripts/pars_hit_qc.py "
        "--log {log} "
        "--datatype {params.datatype} "
        "--timestamp {params.timestamp} "
        "--channel {params.channel} "
        "--configs {configs} "
        "--plot_path {output.plot_file} "
        "--save_path {output.qc_file} "
        "--pulser_file {input.pulser} "
        "--cal_files {input.files} "
        "--fft_files {input.fft_files} "


# This rule builds the energy calibration using the calibration dsp files
rule build_energy_calibration:
    input:
        files=os.path.join(
            filelist_path(setup), "all-{experiment}-{period}-{run}-cal-dsp.filelist"
        ),
        pulser=get_pattern_pars_tmp_channel(setup, "tcm", "pulser_ids"),
        ctc_dict=ancient(
            lambda wildcards: pars_catalog.get_par_file(
                setup, wildcards.timestamp, "dsp"
            )
        ),
        inplots=get_pattern_plts_tmp_channel(setup, "hit", "qc"),
        in_hit_dict=get_pattern_pars_tmp_channel(setup, "hit", "qc"),
    params:
        timestamp="{timestamp}",
        datatype="cal",
        channel="{channel}",
    output:
        ecal_file=temp(get_pattern_pars_tmp_channel(setup, "hit", "energy_cal")),
        results_file=temp(
            get_pattern_pars_tmp_channel(
                setup, "hit", "energy_cal_objects", extension="pkl"
            )
        ),
        plot_file=temp(get_pattern_plts_tmp_channel(setup, "hit", "energy_cal")),
    log:
        get_pattern_log_channel(setup, "pars_hit_energy_cal"),
    group:
        "par-hit"
    resources:
        runtime=300,
    shell:
        "{swenv} python3 -B "
        "{basedir}/../scripts/pars_hit_ecal.py "
        "--log {log} "
        "--datatype {params.datatype} "
        "--timestamp {params.timestamp} "
        "--channel {params.channel} "
        "--configs {configs} "
        "--metadata {meta} "
        "--plot_path {output.plot_file} "
        "--results_path {output.results_file} "
        "--save_path {output.ecal_file} "
        "--inplot_dict {input.inplots} "
        "--in_hit_dict {input.in_hit_dict} "
        "--ctc_dict {input.ctc_dict} "
        "--pulser_file {input.pulser} "
        "--files {input.files}"


# This rule builds the a/e calibration using the calibration dsp files
rule build_aoe_calibration:
    input:
        files=os.path.join(
            filelist_path(setup), "all-{experiment}-{period}-{run}-cal-dsp.filelist"
        ),
        pulser=get_pattern_pars_tmp_channel(setup, "tcm", "pulser_ids"),
        ecal_file=get_pattern_pars_tmp_channel(setup, "hit", "energy_cal"),
        eres_file=get_pattern_pars_tmp_channel(
            setup, "hit", "energy_cal_objects", extension="pkl"
        ),
        inplots=get_pattern_plts_tmp_channel(setup, "hit", "energy_cal"),
    params:
        timestamp="{timestamp}",
        datatype="cal",
        channel="{channel}",
    output:
        hit_pars=temp(get_pattern_pars_tmp_channel(setup, "hit", "aoe_cal")),
        aoe_results=temp(
            get_pattern_pars_tmp_channel(
                setup, "hit", "aoe_cal_objects", extension="pkl"
            )
        ),
        plot_file=temp(get_pattern_plts_tmp_channel(setup, "hit", "aoe_cal")),
    log:
        get_pattern_log_channel(setup, "pars_hit_aoe_cal"),
    group:
        "par-hit"
    resources:
        runtime=300,
    shell:
        "{swenv} python3 -B "
        "{basedir}/../scripts/pars_hit_aoe.py "
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
        "--pulser_file {input.pulser} "
        "--ecal_file {input.ecal_file} "
        "{input.files}"


# This rule builds the lq calibration using the calibration dsp files
rule build_lq_calibration:
    input:
        files=os.path.join(
            filelist_path(setup), "all-{experiment}-{period}-{run}-cal-dsp.filelist"
        ),
        pulser=get_pattern_pars_tmp_channel(setup, "tcm", "pulser_ids"),
        ecal_file=get_pattern_pars_tmp_channel(setup, "hit", "aoe_cal"),
        eres_file=get_pattern_pars_tmp_channel(
            setup, "hit", "aoe_cal_objects", extension="pkl"
        ),
        inplots=get_pattern_plts_tmp_channel(setup, "hit", "aoe_cal"),
    params:
        timestamp="{timestamp}",
        datatype="cal",
        channel="{channel}",
    output:
        hit_pars=temp(get_pattern_pars_tmp_channel(setup, "hit")),
        lq_results=temp(
            get_pattern_pars_tmp_channel(setup, "hit", "objects", extension="pkl")
        ),
        plot_file=temp(get_pattern_plts_tmp_channel(setup, "hit")),
    log:
        get_pattern_log_channel(setup, "pars_hit_lq_cal"),
    group:
        "par-hit"
    resources:
        runtime=300,
    shell:
        "{swenv} python3 -B "
        "{basedir}/../scripts/pars_hit_lq.py "
        "--log {log} "
        "--configs {configs} "
        "--datatype {params.datatype} "
        "--timestamp {params.timestamp} "
        "--inplots {input.inplots} "
        "--channel {params.channel} "
        "--lq_results {output.lq_results} "
        "--eres_file {input.eres_file} "
        "--hit_pars {output.hit_pars} "
        "--plot_file {output.plot_file} "
        "--pulser_file {input.pulser} "
        "--ecal_file {input.ecal_file} "
        "{input.files}"


rule build_pars_hit_objects:
    input:
        lambda wildcards: get_par_chanlist(
            setup,
            f"all-{wildcards.experiment}-{wildcards.period}-{wildcards.run}-cal-{wildcards.timestamp}-channels",
            "hit",
            basedir,
            configs,
            chan_maps,
            name="objects",
            extension="pkl",
        ),
    output:
        get_pattern_pars(
            setup,
            "hit",
            name="objects",
            extension="dir",
            check_in_cycle=check_in_cycle,
        ),
    params:
        ro_input=lambda _, input: ro(input),
    group:
        "merge-hit"
    shell:
        "{swenv} python3 -B "
        "{basedir}/../scripts/merge_channels.py "
        "--input {params.ro_input} "
        "--output {output} "


rule build_plts_hit:
    input:
        lambda wildcards: get_plt_chanlist(
            setup,
            f"all-{wildcards.experiment}-{wildcards.period}-{wildcards.run}-cal-{wildcards.timestamp}-channels",
            "hit",
            basedir,
            configs,
            chan_maps,
        ),
    output:
        get_pattern_plts(setup, "hit"),
    params:
        ro_input=lambda _, input: ro(input),
    group:
        "merge-hit"
    shell:
        "{swenv} python3 -B "
        "{basedir}/../scripts/merge_channels.py "
        "--input {params.ro_input} "
        "--output {output} "


rule build_pars_hit:
    input:
        infiles=lambda wildcards: get_par_chanlist(
            setup,
            f"all-{wildcards.experiment}-{wildcards.period}-{wildcards.run}-cal-{wildcards.timestamp}-channels",
            "hit",
            basedir,
            configs,
            chan_maps,
        ),
        plts=get_pattern_plts(setup, "hit"),
        objects=get_pattern_pars(
            setup,
            "hit",
            name="objects",
            extension="dir",
            check_in_cycle=check_in_cycle,
        ),
    params:
        ro_input=lambda _, input: {k: ro(v) for k, v in input.items()},
    output:
        get_pattern_pars(setup, "hit", check_in_cycle=check_in_cycle),
    group:
        "merge-hit"
    shell:
        "{swenv} python3 -B "
        "{basedir}/../scripts/merge_channels.py "
        "--input {params.ro_input[infiles]} "
        "--output {output} "


rule build_hit:
    input:
        dsp_file=get_pattern_tier_dsp(setup),
        pars_file=lambda wildcards: pars_catalog.get_par_file(
            setup, wildcards.timestamp, "hit"
        ),
    output:
        tier_file=get_pattern_tier(setup, "hit", check_in_cycle=check_in_cycle),
        db_file=get_pattern_pars_tmp(setup, "hit_db"),
    params:
        timestamp="{timestamp}",
        datatype="{datatype}",
        tier="hit",
        ro_input=lambda _, input: {k: ro(v) for k, v in input.items()},
    log:
        get_pattern_log(setup, "tier_hit"),
    group:
        "tier-hit"
    resources:
        runtime=300,
    shell:
        "{swenv} python3 -B "
        "{basedir}/../scripts/build_hit.py "
        f"--configs {ro(configs)} "
        "--log {log} "
        "--tier {params.tier} "
        "--datatype {params.datatype} "
        "--timestamp {params.timestamp} "
        "--pars_file {params.ro_input[pars_file]} "
        "--output {output.tier_file} "
        "--input {params.ro_input[dsp_file]} "
        "--db_file {output.db_file}"
