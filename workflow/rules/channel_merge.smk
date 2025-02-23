import inspect

from legenddataflow import patterns
from legenddataflow.utils import set_last_rule_name
from legenddataflow.execenv import execenv_pyexe


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
        output:
            patterns.get_pattern_plts(config, tier),
        group:
            f"merge-{tier}"
        shell:
            execenv_pyexe(config, "merge-channels") + \
            "--input {input} "
            "--output {output} "

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
        output:
            patterns.get_pattern_pars(
                config,
                tier,
                name="objects",
                extension="dir",
                check_in_cycle=check_in_cycle,
            ),
        group:
            f"merge-{tier}"
        shell:
            execenv_pyexe(config, "merge-channels") + \
            "--input {input} "
            "--output {output} "

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
            output:
                temp(
                    patterns.get_pattern_pars_tmp(
                        config,
                        tier,
                        datatype="cal",
                    )
                ),
            group:
                f"merge-{tier}"
            shell:
                execenv_pyexe(config, "merge-channels") + \
                "--input {input} "
                "--output {output} "

        set_last_rule_name(workflow, f"build_pars_{tier}_db")

        rule:
            """Merge pars for SiPM channels in a single pars file."""
            input:
                lambda wildcards: get_par_chanlist(
                    config,
                    f"all-{wildcards.experiment}-{wildcards.period}-{wildcards.run}-{wildcards.datatype}-{wildcards.timestamp}-channels",
                    tier,
                    basedir,
                    det_status,
                    chan_maps,
                    datatype=wildcards.datatype,
                    system="spms"
                ),
            output:
                patterns.get_pattern_pars(
                    config,
                    tier,
                    datatype="{datatype}",
                ),
            group:
                f"merge-{tier}"
            shell:
                execenv_pyexe(config, "merge-channels") + \
                "--input {input} "
                "--output {output} "

        set_last_rule_name(workflow, f"build_pars_spms_{tier}_db")

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
            in_db=patterns.get_pattern_pars_tmp(
                config,
                tier,
                datatype="cal",
            ) if lh5_merge is True else [],
            plts=patterns.get_pattern_plts(config, tier),
            objects=patterns.get_pattern_pars(
                config,
                tier,
                name="objects",
                extension="dir",
                check_in_cycle=check_in_cycle,
            ),
        output:
            out_file=patterns.get_pattern_pars(
                config,
                tier,
                extension="lh5" if lh5_merge is True else inspect.signature(patterns.get_pattern_pars).parameters['extension'].default,
                check_in_cycle=check_in_cycle,
            ),
            out_db=patterns.get_pattern_pars(config, tier, check_in_cycle=check_in_cycle) if lh5_merge is True else [],
        group:
            f"merge-{tier}"
        run:
            shell_string = (
                execenv_pyexe(config, "merge-channels") + \
                "--output {output.out_file} "
                "--input {input.in_files} "
            )
            if lh5_merge is True:
                shell_string += (
                    "--in-db {input.in_db} "
                    "--out-db {output.out_db} "
                )
            shell(shell_string)

    set_last_rule_name(workflow, f"build_pars_{tier}")
