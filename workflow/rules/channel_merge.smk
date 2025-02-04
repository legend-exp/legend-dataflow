from legenddataflow.patterns import (
    get_pattern_pars_tmp_channel,
    get_pattern_plts_tmp_channel,
    get_pattern_plts,
    get_pattern_tier,
    get_pattern_pars_tmp,
    get_pattern_pars,
)
from legenddataflow.utils import set_last_rule_name
import inspect
from legenddataflow.execenv import execenv_smk_py_script

def build_merge_rules(tier, lh5_merge=False, lh5_tier=None):
    if lh5_tier is None:
        lh5_tier = tier
    rule:
        input:
            lambda wildcards: get_plt_chanlist(
                config,
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
            get_pattern_plts(config, tier),
        group:
            f"merge-{tier}"
        shell:
            f'{execenv_smk_py_script(config, "merge_channels")}'
            "--input {input} "
            "--output {output} "
            "--channelmap {meta} "

    set_last_rule_name(workflow, f"build_plts_{tier}")

    rule:
        input:
            lambda wildcards: get_par_chanlist(
                config,
                f"all-{wildcards.experiment}-{wildcards.period}-{wildcards.run}-cal-{wildcards.timestamp}-channels",
                tier,
                basedir,
                det_status,
                chan_maps,
                name="objects",
                extension="pkl",
            ),
        params:
            timestamp="{timestamp}",
            datatype="cal",
        output:
            get_pattern_pars(
                config,
                tier,
                name="objects",
                extension="dir",
                check_in_cycle=check_in_cycle,
            ),
        group:
            f"merge-{tier}"
        shell:
            f'{execenv_smk_py_script(config, "merge_channels")}'
            "--input {input} "
            "--output {output} "
            "--timestamp {params.timestamp} "
            "--channelmap {meta} "

    set_last_rule_name(workflow, f"build_pars_{tier}_objects")

    if lh5_merge is True:
        rule:
            input:
                lambda wildcards: get_par_chanlist(
                    config,
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
                        config,
                        tier,
                        datatype="cal",
                    )
                ),
            group:
                f"merge-{tier}"
            shell:
                f'{execenv_smk_py_script(config, "merge_channels")}'
                "--input {input} "
                "--output {output} "
                "--timestamp {params.timestamp} "
                "--channelmap {meta} "

        set_last_rule_name(workflow, f"build_pars_{tier}_db")

    rule:
        input:
            in_files=lambda wildcards: get_par_chanlist(
                config,
                f"all-{wildcards.experiment}-{wildcards.period}-{wildcards.run}-cal-{wildcards.timestamp}-channels",
                lh5_tier,
                basedir,
                det_status,
                chan_maps,
                extension="lh5" if lh5_merge is True else inspect.signature(get_par_chanlist).parameters['extension'].default,
            ),
            in_db=get_pattern_pars_tmp(
                config,
                tier,
                datatype="cal",
            ) if lh5_merge is True else [],
            plts=get_pattern_plts(config, tier),
            objects=get_pattern_pars(
                config,
                tier,
                name="objects",
                extension="dir",
                check_in_cycle=check_in_cycle,
            ),
        params:
            timestamp="{timestamp}",
            datatype="cal",
        output:
            out_file=get_pattern_pars(
                config,
                tier,
                extension="lh5" if lh5_merge is True else inspect.signature(get_pattern_pars).parameters['extension'].default,
                check_in_cycle=check_in_cycle,
            ),
            out_db=get_pattern_pars(config, tier, check_in_cycle=check_in_cycle) if lh5_merge is True else [],
        group:
            f"merge-{tier}"
        run:
            shell_string = (
                f'{execenv_smk_py_script(config, "merge_channels")}'
                "--output {output.out_file} "
                "--input {input.in_files} "
                "--timestamp {params.timestamp} "
                "--channelmap {meta} "
            )
            if lh5_merge is True:
                shell_string += (
                    "--in_db {input.in_db} "
                    "--out_db {output.out_db} "
                )
            shell(shell_string)

    set_last_rule_name(workflow, f"build_pars_{tier}")
