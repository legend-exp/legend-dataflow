import glob
import json
import os

from scripts.util.FileKey import FileKey, run_grouper
from scripts.util.patterns import get_pattern_tier, get_pattern_tier_raw_blind


def get_analysis_runs(ignore_keys_file=None, analysis_runs_file=None):
    ignore_keys = []
    if ignore_keys_file is not None:
        if os.path.isfile(ignore_keys_file):
            with open(ignore_keys_file) as f:
                ignore_keys = f.read().splitlines()
            ignore_keys = [
                key.split("#")[0].strip() if "#" in key else key.strip()
                for key in ignore_keys
            ]
        else:
            print("no ignore_keys.keylist file found")
            ignore_keys = []
    else:
        ignore_keys = []

    if analysis_runs_file is not None:
        if os.path.isfile(analysis_runs_file):
            with open(analysis_runs_file) as f:
                analysis_runs = json.load(f)
        else:
            analysis_runs = []
            print("no analysis_runs file found")
    else:
        analysis_runs = []
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


def get_pattern(setup, tier):
    if tier == "blind":
        fn_pattern = get_pattern_tier(setup, "raw", check_in_cycle=False)
    elif tier == "skm" or tier == "pet_concat":
        fn_pattern = get_pattern_tier(setup, "pet", check_in_cycle=False)
    elif tier == "evt_concat":
        fn_pattern = get_pattern_tier(setup, "evt", check_in_cycle=False)
    else:
        fn_pattern = get_pattern_tier(setup, tier, check_in_cycle=False)
    return fn_pattern


def concat_phy_filenames(setup, phy_filenames, tier):
    fn_pattern = get_pattern(setup, tier)
    sorted_phy_filenames = run_grouper(phy_filenames)
    phy_filenames = []
    for run in sorted_phy_filenames:
        key = FileKey.get_filekey_from_pattern(run[0], fn_pattern)
        out_key = FileKey.get_path_from_filekey(
            key, get_pattern_tier(setup, tier, check_in_cycle=False)
        )[0]

        phy_filenames.append(out_key)
    return phy_filenames


def build_filelist(
    setup,
    filekeys,
    search_pattern,
    tier,
    ignore_keys=None,
    analysis_runs=None,
    file_selection="all",
):
    fn_pattern = get_pattern(setup, tier)

    if ignore_keys is None:
        ignore_keys = []
    if analysis_runs is None:
        analysis_runs = []

    phy_filenames = []
    other_filenames = []

    for key in filekeys:
        fn_glob_pattern = key.get_path_from_filekey(search_pattern)[0]
        files = glob.glob(fn_glob_pattern)
        for f in files:
            _key = FileKey.get_filekey_from_pattern(f, search_pattern)
            if _key.name in ignore_keys:
                pass
            else:
                if tier == "blind" and _key.datatype == "phy":
                    filename = FileKey.get_path_from_filekey(
                        _key, get_pattern_tier_raw_blind(setup)
                    )
                elif tier == "skm":
                    filename = FileKey.get_path_from_filekey(
                        _key, get_pattern_tier(setup, "pet", check_in_cycle=False)
                    )
                else:
                    filename = FileKey.get_path_from_filekey(_key, fn_pattern)

                if file_selection == "all":
                    if _key.datatype == "phy":
                        phy_filenames += filename
                    else:
                        other_filenames += filename
                elif file_selection == "sel":
                    if analysis_runs == "all" or (
                        _key.period in analysis_runs
                        and (
                            _key.run in analysis_runs[_key.period]
                            or analysis_runs[_key.period] == "all"
                        )
                    ):
                        if _key.datatype == "phy":
                            phy_filenames += filename
                        else:
                            other_filenames += filename
                else:
                    msg = "unknown file selection"
                    raise ValueError(msg)

    phy_filenames = sorted(phy_filenames)
    other_filenames = sorted(other_filenames)
    if tier == "skm" or tier == "pet_concat" or tier == "evt_concat":
        concat_phy_filenames(setup, phy_filenames, tier)
    return phy_filenames + other_filenames


def get_filelist(
    wildcards, setup, search_pattern, ignore_keys_file=None, analysis_runs_file=None
):
    file_selection = wildcards.label[:3]
    keypart = wildcards.label[3:]

    analysis_runs, ignore_keys = get_analysis_runs(ignore_keys_file, analysis_runs_file)

    filekeys = get_keys(keypart)
    return build_filelist(
        setup,
        filekeys,
        search_pattern,
        wildcards.tier,
        ignore_keys,
        analysis_runs,
        file_selection,
    )


def get_filelist_full_wildcards(
    wildcards,
    setup,
    search_pattern,
    tier,
    ignore_keys_file=None,
    analysis_runs_file=None,
    file_selection="all",
):
    keypart = f"-{wildcards.experiment}-{wildcards.period}-{wildcards.run}-{wildcards.datatype}"

    analysis_runs, ignore_keys = get_analysis_runs(ignore_keys_file, analysis_runs_file)

    filekeys = get_keys(keypart)
    return build_filelist(
        setup,
        filekeys,
        search_pattern,
        tier,
        ignore_keys,
        analysis_runs,
        file_selection,
    )
