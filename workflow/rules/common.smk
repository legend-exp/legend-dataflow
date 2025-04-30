"""
Helper functions for running data production
"""

from pathlib import Path
import json
import re
from legenddataflow import patterns as patt
from dbetto.catalog import Catalog
from dbetto import TextDB
from legenddataflow import utils


def ro(path):
    return utils.as_ro(config, path)


def get_blinding_curve_file(wildcards):
    """func to get the blinding calibration curves from the overrides"""
    par_files = Catalog.get_files(
        Path(patt.par_overwrite_path(config)) / "raw" / "validity.yaml",
        wildcards.timestamp,
    )
    if isinstance(par_files, str):
        return str(Path(patt.par_overwrite_path(config)) / "raw" / par_files)
    else:
        return [
            str(Path(patt.par_overwrite_path(config)) / "raw" / par_file)
            for par_file in par_files
        ]


def get_blinding_check_file(
    wildcards,
):
    """func to get the right blinding check file"""
    par_files = raw_catalog.get_files(wildcards.timestamp)
    if isinstance(par_files, str):
        return Path(patt.get_pars_path(config, "raw")) / par_files
    else:
        return [
            Path(patt.get_pars_path(config, "raw")) / par_file for par_file in par_files
        ]


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


def get_input_par_file(config, wildcards, tier, name):
    allow_none = config.get("allow_none", False)
    par_overwrite_file = Path(patt.par_overwrite_path(config)) / tier / "validity.yaml"
    pars_files_overwrite = Catalog.get_files(
        par_overwrite_file,
        wildcards.timestamp,
        category=wildcards.datatype if hasattr(wildcards, "datatype") else "all",
    )
    for pars_file in pars_files_overwrite:
        if name in str(pars_file):
            return Path(patt.par_overwrite_path(config)) / tier / pars_file
    if allow_none or (hasattr(wildcards, "datatype") & (wildcards.datatype != "phy")):
        return []
    else:
        raise ValueError(f"Could not find model in {pars_files_overwrite}")


def get_overwrite_file(config, tier, wildcards=None, timestamp=None, name=None):
    par_overwrite_file = Path(patt.par_overwrite_path(config)) / tier / "validity.yaml"
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
            out_files.append(Path(patt.par_overwrite_path(config)) / tier / pars_file)
    if len(out_files) == 0:
        raise ValueError(f"Could not find name in {pars_files_overwrite}")
    else:
        return out_files


def get_search_pattern(config, tier):
    """
    This func gets the search pattern for the relevant tier passed.
    """
    if tier == "daq":
        return patt.get_pattern_tier_daq_unsorted(config, extension="*")
    elif tier == "raw":
        return patt.get_pattern_tier_daq(config, extension="*")
    else:
        return patt.get_pattern_tier(config, "raw", check_in_cycle=False)


def get_table_name(metadata, config, datatype, timestamp, detector, tier):
    if isinstance(metadata, (str, Path)):
        chmap = metadata.channelmap(timestamp, system=datatype)
    elif isinstance(metadata, Catalog):
        chmap = metadata.valid_for(timestamp, system=datatype)
    else:
        raise ValueError(
            f"metadata must be a string or a Catalog object, not {type(metadata)}"
        )
    return config.table_format[tier].format(ch=chmap[detector].daq.rawid)


def get_all_channels(channelmap, timestamp, datatype):
    if isinstance(channelmap, (str, Path)):
        channelmap = TextDB(channelmap, lazy=True)

    if isinstance(channelmap, TextDB):
        chmap = channelmap.channelmaps.on(timestamp, system=datatype)
    else:
        chmap = channelmap.valid_for(timestamp, system=datatype)

    channels = list(chmap)

    if len(channels) == 0:
        print("WARNING: No channels found")  # noqa: T201

    return channels


def get_table_mapping(channelmap, timestamp, datatype, tier):
    if isinstance(channelmap, (str, Path)):
        channelmap = TextDB(channelmap, lazy=True)
    channel_dict = channelmap.valid_for(timestamp, system=datatype)
    detectors = get_all_channels(channelmap, timestamp, datatype)
    return json.dumps(
        {
            detector: f"ch{channel_dict[detector].daq.rawid:07}/{tier}"
            for detector in detectors
        }
    )


def get_alias(channelmap, timestamp, datatype, tier):
    if isinstance(channelmap, (str, Path)):
        channelmap = TextDB(channelmap, lazy=True)
    channel_dict = channelmap.valid_for(timestamp, system=datatype)
    detectors = get_all_channels(channelmap, timestamp, datatype)
    return json.dumps(
        {
            f"ch{channel_dict[detector].daq.rawid:07}/{tier}": f"{detector}/{tier}"
            for detector in detectors
        }
    )


def strip_channel_wildcard_constraint(files):
    if isinstance(files, str):
        files = [files]
    return [re.sub(r"\{channel,[^}]+\}", "{channel}", file) for file in files]


def get_threads(wildcards):
    if config.get("multiprocess", False):
        mode = config.get("multiprocess_mode", "max_usage")
        if mode == "max_usage":
            return min(config.get("max_processes", 50), workflow.cores)
        elif mode == "single":
            return min(
                len(
                    get_table_mapping(
                        channelmap_textdb,
                        wildcards.timestamp,
                        wildcards.datatype,
                        "raw",
                    )
                ),
                workflow.cores,
                64,
            )
    else:
        return 1


def _make_input_pars_file(wildcards):
    """Prepare the input pars files for the `build_dsp` rule."""
    # first get the files from the catalog
    filelist = dsp_par_catalog.get_par_file(config, wildcards.timestamp, "dsp")

    # then add the spms par files
    if wildcards.datatype not in ("cal", "xtc"):
        filelist += [
            patt.get_pattern_pars(config, "dsp", name="spms", datatype="{datatype}")
        ]

    return filelist
