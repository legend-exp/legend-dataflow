from scripts.util.patterns import (
    get_pattern_pars_tmp_channel,
    get_pattern_plts_tmp_channel,
    get_pattern_plts,
    get_pattern_tier,
    get_pattern_pars_tmp,
    get_pattern_pars,
)
from scripts.util.utils import set_last_rule_name
import inspect

def build_merge_rules(tier,lh5_merge=False):
    rule:
        input:
            lambda wildcards: get_plt_chanlist(
                setup,
                f"all-{wildcards.experiment}-{wildcards.period}-{wildcards.run}-cal-{wildcards.timestamp}-channels",
                tier,
                basedir,
                det_status,
                chan_maps,
            ),
        params:
            timestamp="{timestamp}",
            datatype="cal",
        output:
            get_pattern_plts(setup, tier),
        group:
            f"merge-{tier}"
        shell:
            "{swenv} python3 -B "
            "{basedir}/../scripts/merge_channels.py "
            "--input {input} "
            "--output {output} "
            "--channelmap {meta} "

    set_last_rule_name(workflow, f"build_plts_{tier}")

    rule:
        input:
            lambda wildcards: get_par_chanlist(
                setup,
                f"all-{wildcards.experiment}-{wildcards.period}-{wildcards.run}-cal-{wildcards.timestamp}-channels",
                tier,
                basedir,
                det_status,
                chan_maps,
                name="objects",
                extension="pkl",
            ),
        output:
            get_pattern_pars(
                setup,
                tier,
                name="objects",
                extension="dir",
                check_in_cycle=check_in_cycle,
            ),
        group:
            f"merge-{tier}"
        shell:
            "{swenv} python3 -B "
            "{basedir}/../scripts/merge_channels.py "
            "--input {input} "
            "--output {output} "
            "--timestamp {params.timestamp} "
            "--channelmap {meta} "

    set_last_rule_name(workflow, f"build_pars_{tier}_objects")

    if lh5_merge is True:
        rule:
            input:
                lambda wildcards: get_par_chanlist(
                    setup,
                    f"all-{wildcards.experiment}-{wildcards.period}-{wildcards.run}-cal-{wildcards.timestamp}-channels",
                    tier,
                    basedir,
                    det_status,
                    chan_maps,
                ),
            params:
                timestamp="{timestamp}",
                datatype="cal",
            output:
                temp(
                    get_pattern_pars_tmp(
                        setup,
                        tier,
                        datatype="cal",
                    )
                ),
            group:
                f"merge-{tier}"
            shell:
                "{swenv} python3 -B "
                "{basedir}/../scripts/merge_channels.py "
                "--input {input} "
                "--output {output} "
                "--timestamp {params.timestamp} "
                "--channelmap {meta} "

        set_last_rule_name(workflow, f"build_pars_{tier}_db")

    rule:
        input:
            in_files=lambda wildcards: get_par_chanlist(
                setup,
                f"all-{wildcards.experiment}-{wildcards.period}-{wildcards.run}-cal-{wildcards.timestamp}-channels",
                tier,
                basedir,
                det_status,
                chan_maps,
                extension="lh5" if lh5_merge is True else inspect.signature(get_par_chanlist).parameters['extension'].default,
            ),
            in_db=get_pattern_pars_tmp(
                setup,
                "dsp",
                datatype="cal",
            ) if lh5_merge is True else None,
            plts=get_pattern_plts(setup, "dsp"),
            objects=get_pattern_pars(
                setup,
                "dsp",
                name="objects",
                extension="dir",
                check_in_cycle=check_in_cycle,
            ),
        params:
            timestamp="{timestamp}",
            datatype="cal",
        output:
            out_file=get_pattern_pars(
                setup,
                tier,
                extension="lh5" if lh5_merge is True else inspect.signature(get_pattern_pars).parameters['extension'].default,
                check_in_cycle=check_in_cycle,
            ),
            out_db=get_pattern_pars(setup, tier, check_in_cycle=check_in_cycle) if lh5_merge is True else None,
        group:
            f"merge-{tier}"
        run:
            shell_cmd  = "{swenv} python3 -B "
            shell_cmd += "{basedir}/../scripts/merge_channels.py "
            shell_cmd += "--output {output.out_file} "
            shell_cmd += "--input {input.in_files} "
            shell_cmd += "--timestamp {params.timestamp} "
            shell_cmd += "--channelmap {meta} "
            if lh5_merge is True:
                shell_cmd +="--in_db {input.in_db} "
                shell_cmd +="--out_db {output.out_db} "
            shell(
                shell_cmd
            )

    set_last_rule_name(workflow, f"build_pars_{tier}")
