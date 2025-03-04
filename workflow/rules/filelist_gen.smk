import glob
import json, yaml
from pathlib import Path

from legenddataflow.FileKey import FileKey, run_grouper
from legenddataflow import patterns as patt

concat_datatypes = ["phy"]
concat_tiers = ["skm", "pet_concat", "evt_concat"]
blind_datatypes = ["phy"]


def expand_runs(in_dict):
    """
    This function expands out the runs if a range is specified in the dictionary
    e.g.
    {
        "p01": "r001..r005"
    }
    """
    for per, datalist in in_dict.items():
        for datatype, run_list in datalist.items():
            if isinstance(run_list, str) and ".." in runs:
                start, end = runs.split("..")
                in_dict[per][datatype] = [
                    f"r{x:03}" for x in range(int(start[1:]), int(end[1:]) + 1)
                ]
    return in_dict


def get_analysis_runs(
    ignore_keys_file=None, analysis_runs_file=None, file_selection="all"
):
    """
    This function reads in the ignore_keys and analysis_runs files and returns the dictionaries
    """
    ignore_keys = []
    analysis_runs = {}
    if ignore_keys_file is not None:
        if Path(ignore_keys_file).is_file():
            if Path(ignore_keys_file).suffix == ".json":
                with Path(ignore_keys_file).open() as f:
                    ignore_keys = json.load(f)
            elif Path(ignore_keys_file).suffix == ".keylist":
                with Path(ignore_keys_file).open() as f:
                    ignore_keys = f.read().splitlines()
            elif Path(ignore_keys_file).suffix in (".yaml", ".yml"):
                with Path(ignore_keys_file).open() as f:
                    ignore_keys = yaml.safe_load(f)
            else:
                raise ValueError(
                    "ignore_keys_file file not in json, yaml or keylist format"
                )
            ignore_keys = [  # remove any comments in the keylist
                key.split("#")[0].strip() if "#" in key else key.strip()
                for key in ignore_keys
            ]
        else:
            msg = f"no ignore_keys file found: {ignore_keys_file}"
            raise ValueError(msg)

    if analysis_runs_file is not None and file_selection != "all":
        if Path(analysis_runs_file).is_file():
            if Path(ignore_keys_file).suffix == ".json":
                with Path(analysis_runs_file).open() as f:
                    analysis_runs = json.load(f)
            elif Path(ignore_keys_file).suffix in (".yaml", ".yml"):
                with Path(analysis_runs_file).open() as f:
                    analysis_runs = yaml.safe_load(f)
            else:
                msg = f"analysis_runs file not in json or yaml format: {analysis_runs_file}"
                raise ValueError(msg)
            if file_selection in analysis_runs:
                analysis_runs = expand_runs(
                    analysis_runs[file_selection]
                )  # select the file_selection and expand out the runs
            else:
                msg = f"Unknown file selection: {file_selection} not in {list(analysis_runs)}"
                raise ValueError(msg)
        else:
            msg = f"no analysis_runs file found: {analysis_runs_file}"
            raise ValueError(msg)
    return analysis_runs, ignore_keys


def get_keys(keypart):

    key = FileKey.parse_keypart(keypart)

    item_list = []
    for item in key:
        _item = item.split("_") if "_" in item else item
        if isinstance(_item, list):
            item_list.append(_item)
        else:
            item_list.append([_item])

    filekeys = []
    for i in item_list[0]:
        for j in item_list[1]:
            for k in item_list[2]:
                for i2 in item_list[3]:
                    for j2 in item_list[4]:
                        filekeys.append(FileKey(i, j, k, i2, j2))
    return filekeys


def get_pattern(config, tier):
    """
    Helper function to get the search pattern for the given tier,
    some tiers such as skm need to refer to a different pattern when looking for files
    as only phy files are taken to skm others are only taken to pet
    """
    if tier == "blind":
        fn_pattern = patt.get_pattern_tier(config, "raw", check_in_cycle=False)
    elif tier in ("skm", "pet_concat"):
        fn_pattern = patt.get_pattern_tier(config, "pet", check_in_cycle=False)
    elif tier == "evt_concat":
        fn_pattern = patt.get_pattern_tier(config, "evt", check_in_cycle=False)
    elif tier == "daq":
        fn_pattern = patt.get_pattern_tier_daq(config, extension="{ext}")
    else:
        fn_pattern = patt.get_pattern_tier(config, tier, check_in_cycle=False)
    return fn_pattern


def concat_phy_filenames(config, phy_filenames, tier):
    """
    This function concatenates the files from the same run together
    """
    fn_pattern = get_pattern(config, tier)
    # group files by run
    sorted_phy_filenames = run_grouper(phy_filenames)
    phy_filenames = []

    for run in sorted_phy_filenames:
        key = FileKey.get_filekey_from_pattern(run[0], fn_pattern)
        out_key = FileKey.get_path_from_filekey(
            key, patt.get_pattern_tier(config, tier, check_in_cycle=False)
        )[0]

        phy_filenames.append(out_key)

    return phy_filenames


def build_filelist(
    config,
    filekeys,
    search_pattern,
    tier,
    ignore_keys=None,
    analysis_runs=None,
):
    """
    This function builds the filelist for the given filekeys, search pattern
    and tier. It will ignore any keys in the ignore_keys list and only include
    the keys specified in the analysis_runs dict.
    """
    fn_pattern = get_pattern(config, tier)

    if ignore_keys is None:
        ignore_keys = []
    if analysis_runs is None:
        analysis_runs = {}

    phy_filenames = []
    other_filenames = []

    for key in filekeys:
        if Path(search_pattern).suffix == ".*":
            search_pattern = Path(search_pattern).with_suffix(".{ext}")
        fn_glob_pattern = key.get_path_from_filekey(search_pattern, ext="*")[0]
        files = glob.glob(fn_glob_pattern)
        for f in files:
            _key = FileKey.get_filekey_from_pattern(f, search_pattern)
            if _key.name in ignore_keys:
                pass
            else:
                if tier == "blind" and _key.datatype in blind_datatypes:
                    filename = FileKey.get_path_from_filekey(
                        _key, patt.get_pattern_tier_raw_blind(config)
                    )
                elif tier == "skm":
                    filename = FileKey.get_path_from_filekey(
                        _key, patt.get_pattern_tier(config, "pet", check_in_cycle=False)
                    )
                elif tier == "daq":
                    filename = FileKey.get_path_from_filekey(
                        _key, fn_pattern.with_suffix(Path(f).suffix)
                    )
                else:
                    filename = FileKey.get_path_from_filekey(_key, fn_pattern)

                if analysis_runs == {}:
                    if (
                        _key.datatype in concat_datatypes
                    ):  # separate out phy files as some tiers these are concatenated
                        phy_filenames += filename
                    else:
                        other_filenames += filename
                else:
                    # check if period in analysis_runs dicts
                    # check if run in analysis_runs dicts
                    # or if runs is just specified as "all"
                    if (
                        _key.datatype in analysis_runs
                        and _key.period in analysis_runs[_key.datatype]
                        and (
                            _key.run in analysis_runs[_key.datatype][_key.period]
                            or analysis_runs[_key.datatype][_key.period] == "all"
                        )
                    ):
                        if _key.datatype in concat_datatypes:
                            phy_filenames += filename  # separate out phy files as some tiers these are concatenated
                        else:
                            other_filenames += filename

    phy_filenames = sorted(phy_filenames)
    other_filenames = sorted(other_filenames)

    if tier in concat_tiers:
        phy_filenames = concat_phy_filenames(
            config, phy_filenames, tier
        )  # concat phy files

    return phy_filenames + other_filenames


def get_filelist(
    wildcards, config, search_pattern, ignore_keys_file=None, analysis_runs_file=None
):
    file_selection = wildcards.label.split("-", 1)[0]
    # remove the file selection from the keypart
    keypart = f'-{wildcards.label.split("-",1)[1]}'
    analysis_runs, ignore_keys = get_analysis_runs(
        ignore_keys_file, analysis_runs_file, file_selection
    )

    filekeys = get_keys(keypart)

    return build_filelist(
        config,
        filekeys,
        search_pattern,
        wildcards.tier,
        ignore_keys,
        analysis_runs,
    )


def get_filelist_full_wildcards(
    wildcards,
    config,
    search_pattern,
    tier,
    ignore_keys_file=None,
    analysis_runs_file=None,
    file_selection="all",
):
    keypart = f"-{wildcards.experiment}-{wildcards.period}-{wildcards.run}-{wildcards.datatype}"

    analysis_runs, ignore_keys = get_analysis_runs(
        ignore_keys_file, analysis_runs_file, file_selection
    )

    filekeys = get_keys(keypart)
    return build_filelist(
        config,
        filekeys,
        search_pattern,
        tier,
        ignore_keys,
        analysis_runs,
    )
