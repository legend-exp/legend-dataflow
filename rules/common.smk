"""
Helper functions for running data production
"""

from pathlib import Path
from scripts.library.patterns import (
    get_pattern_tier_daq_unsorted,
    get_pattern_tier_daq,
    get_pattern_tier,
    par_overwrite_path,
    get_pars_path,
)
from scripts.library import ProcessingFileKey
from dbetto.catalog import Catalog
from scripts.library import utils


def ro(path):
    return utils.as_ro(setup, path)


def get_blinding_curve_file(wildcards):
    """func to get the blinding calibration curves from the overrides"""
    par_files = Catalog.get_files(
        Path(par_overwrite_path(setup)) / "raw" / "validity.yaml",
        wildcards.timestamp,
    )
    if isinstance(par_files, str):
        return str(Path(par_overwrite_path(setup)) / "raw" / par_files)
    else:
        return [
            str(Path(par_overwrite_path(setup)) / "raw" / par_file)
            for par_file in par_files
        ]


def get_blinding_check_file(wildcards):
    """func to get the right blinding check file"""
    par_files = Catalog.get_files(
        Path(get_pars_path(setup, "raw")) / "validity.yaml", wildcards.timestamp
    )
    if isinstance(par_files, str):
        return Path(get_pars_path(setup, "raw")) / par_files
    else:
        return [Path(get_pars_path(setup, "raw")) / par_file for par_file in par_files]


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


def get_input_par_file(wildcards, tier, name):
    par_overwrite_file = Path(par_overwrite_path(setup)) / tier / "validity.yaml"
    pars_files_overwrite = Catalog.get_files(
        par_overwrite_file,
        wildcards.timestamp,
    )
    for pars_file in pars_files_overwrite:
        if name in str(pars_file):
            return Path(par_overwrite_path(setup)) / tier / pars_file
    raise ValueError(f"Could not find model in {pars_files_overwrite}")


def get_overwrite_file(tier, wildcards=None, timestamp=None, name=None):
    par_overwrite_file = Path(par_overwrite_path(setup)) / tier / "validity.yaml"
    if timestamp is not None:
        pars_files_overwrite = Catalog.get_files(
            par_overwrite_file,
            timestamp,
        )
    else:
        pars_files_overwrite = Catalog.get_files(
            par_overwrite_file,
            wildcards.timestamp,
        )
    if name is None:
        fullname = f"{tier}-overwrite.yaml"
    else:
        fullname = f"{tier}_{name}-overwrite.yaml"
    out_files = []
    for pars_file in pars_files_overwrite:
        if fullname in str(pars_file):
            out_files.append(Path(par_overwrite_path(setup)) / tier / pars_file)
    if len(out_files) == 0:
        raise ValueError(f"Could not find name in {pars_files_overwrite}")
    else:
        return out_files


def get_search_pattern(tier):
    """
    This func gets the search pattern for the relevant tier passed.
    """
    if tier == "daq":
        return get_pattern_tier_daq_unsorted(setup, extension="*")
    elif tier == "raw":
        return get_pattern_tier_daq(setup, extension="*")
    else:
        return get_pattern_tier(setup, "raw", check_in_cycle=False)
