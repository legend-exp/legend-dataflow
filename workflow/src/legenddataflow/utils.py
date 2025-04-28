"""
This module contains all the utility needed for the data production.
for substituting the pathvar in the config, also the conversion
from timestamp to unix time
"""

import copy
import os
import re
import string
from pathlib import Path


def subst_vars_impl(x, var_values, ignore_missing=False):
    if isinstance(x, str):
        if "$" in x:
            if ignore_missing:
                return string.Template(x).safe_substitute(var_values)
            else:
                return string.Template(x).substitute(var_values)
        else:
            return x
    if isinstance(x, dict):
        for key in x:
            value = x[key]
            new_value = subst_vars_impl(value, var_values, ignore_missing)
            if new_value is not value:
                x[key] = new_value
        return x
    if isinstance(x, list):
        for i in range(len(x)):
            value = x[i]
            new_value = subst_vars_impl(value, var_values, ignore_missing)
            if new_value is not value:
                x[i] = new_value
        return x
    else:
        return x


def subst_vars(
    props,
    var_values=None,
    use_env=False,
    ignore_missing=False,
):
    if var_values is None:
        var_values = {}
    combined_var_values = var_values
    if use_env:
        combined_var_values = dict(iter(os.environ.items()))
        combined_var_values.update(copy.copy(var_values))
    subst_vars_impl(props, combined_var_values, ignore_missing)


def subst_vars_in_snakemake_config(workflow, config):
    if len(workflow.overwrite_configfiles) == 0:
        msg = "configfile not set!"
        raise RuntimeError(msg)

    config_filename = workflow.overwrite_configfiles[0]

    subst_vars(
        config,
        var_values={"_": Path(config_filename).parent},
        use_env=True,
        ignore_missing=False,
    )
    if "system" in config:
        config["execenv"] = config["execenv"][config["system"]]
    else:
        config["execenv"] = config["execenv"]["bare"]


def set_last_rule_name(workflow, new_name):
    """Sets the name of the most recently created rule to be `new_name`.
    Useful when creating rules dynamically (i.e. unnamed).

    Warning
    -------
    This could mess up the workflow. Use at your own risk.
    """
    rules = workflow._rules
    last_key = next(reversed(rules))
    assert last_key == rules[last_key].name

    rules[new_name] = rules.pop(last_key)
    rules[new_name].name = new_name

    if workflow.default_target == last_key:
        workflow.default_target = new_name

    if last_key in workflow._localrules:
        workflow._localrules.remove(last_key)
        workflow._localrules.add(new_name)

    workflow.check_localrules()


def as_ro(config, path):
    if (
        "read_only_fs_sub_pattern" not in config
        or config["read_only_fs_sub_pattern"] is None
    ):
        return path

    sub_pattern = config["read_only_fs_sub_pattern"]

    if isinstance(path, str):
        return re.sub(*sub_pattern, path)
    if isinstance(path, Path):
        return Path(re.sub(*sub_pattern, path.name))

    return [as_ro(config, p) for p in path]
