import inspect

from legenddataflow.methods import patterns
from legenddataflowscripts.workflow import execenv_pyexe, set_last_rule_name
from legenddataflow.scripts.flow.build_chanlist import get_plt_chanlist, get_par_chanlist


def build_in_channel_merge_rules(rules, tier, output_name=None):


    input_dsp_pars = []
    input_dsp_pars_lh5 = []
    input_dsp_objects = []
    input_dsp_plots = []

    name = f"tier" if output_name is None else f"{tier}_{output_name}"

    group_name = f"merge-inchan-{tier}" 
    if output_name is not None:
        group_name += f"_{output_name}"

    for r in [rules]:
        if hasattr(r.output, dsp_pars):
            input_dsp_pars.append(r.output.dsp_pars)
        if hasattr(r.output, dsp_pars_lh5):
            input_dsp_pars_lh5.append(r.output.dsp_pars_lh5)
        if hasattr(r.output, plots):
            input_dsp_plots.append(r.output.dsp_plots)
        if hasattr(r.output, dsp_objects):
            input_dsp_objects.append(r.output.dsp_objects)

    par_file = temp(get_pattern_pars_tmp_channel(config, "dsp", name=output_name)),
            lh5_file = temp(get_pattern_pars_tmp_channel(config, "dsp", name=output_name, extension="lh5")),
            object_file = 
            plot_file = 

    rule:
        input:
            input_dsp_plots,
        output:
            temp(get_pattern_plts_tmp_channel(config, "dsp", name=output_name)),
        group:
            group_name
        shell:
            execenv_pyexe(config, "merge-in-channels") + \
            "--input {input} "
            "--output {output} "

    set_last_rule_name(workflow, f"build_plts_inchan_{name}")

    rule:
        input:
            input_dsp_objects,
        output:
            temp(get_pattern_pars_tmp_channel(config, "dsp", name=output_name, extension="pkl")),
        group:
            group_name
        shell:
            execenv_pyexe(config, "merge-in-channels") + \
            "--input {input} "
            "--output {output} "

    set_last_rule_name(workflow, f"build_pars_inchan_{name}_objects")

    if len(input_dsp_pars_lh5) > 0:
        rule:
            input:
                input_dsp_pars_lh5,
            output:
                temp(get_pattern_pars_tmp_channel(config, "dsp", name=output_name + "_tmp" if output_name is not None else "tmp")),
            group:
                group_name
            shell:
                execenv_pyexe(config, "merge-in-channels") + \
                "--input {input} "
                "--output {output} "

        set_last_rule_name(workflow, f"build_pars_inchan_{name}_db")


    rule:
        input:
            in_files=input_dsp_pars_lh5 if len(input_dsp_pars_lh5) > 0 else input_dsp_pars,
            in_db=temp(get_pattern_pars_tmp_channel(config, "dsp", name=output_name + "_tmp" if output_name is not None else "tmp")) if len(input_dsp_pars_lh5) > 0 else [],
            plts=temp(get_pattern_plts_tmp_channel(config, "dsp", name=output_name)),
            objects=temp(get_pattern_pars_tmp_channel(config, "dsp", name=output_name, extension="pkl")),
        output:
            out_file=temp(get_pattern_pars_tmp_channel(config, "dsp", name=output_name, extension="lh5")) if len(input_dsp_pars_lh5) > 0 else temp(get_pattern_pars_tmp_channel(config, "dsp", name=output_name)),
            out_db=temp(get_pattern_pars_tmp_channel(config, "dsp", name=output_name)) if len(input_dsp_pars_lh5) > 0 else [],
        group:
            group_name
        run:
            shell_string = (
                execenv_pyexe(config, "merge-in-channels") + \
                "--output {output.out_file} "
                "--input {input.in_files} "
            )
            if len(input_dsp_pars_lh5) > 0:
                shell_string += (
                    "--in-db {input.in_db} "
                    "--out-db {output.out_db} "
                )
            shell(shell_string)

    set_last_rule_name(workflow, f"build_pars_inchan_{name}")
