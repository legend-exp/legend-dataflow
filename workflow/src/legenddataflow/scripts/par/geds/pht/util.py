import pickle as pkl
import re
from pathlib import Path

import numpy as np
from dbetto import Props
from legenddatafow.FileKey import ChannelProcKey, ProcessingFileKey, run_splitter


def update_cal_dicts(cal_dicts, update_dict):
    if re.match(r"(\d{8})T(\d{6})Z", next(iter(cal_dicts))):
        for tstamp in cal_dicts:
            if tstamp in update_dict:
                cal_dicts[tstamp].update(update_dict[tstamp])
            else:
                cal_dicts[tstamp].update(update_dict)
    else:
        cal_dicts.update(update_dict)
    return cal_dicts


def get_run_dict(files):
    out_dicts = {}
    for file in files:
        if Path(file).suffix == ".pkl":
            with Path(file).open("rb") as o:
                vals = pkl.load(o)
        elif Path(file).suffix in (".json", ".yml", ".yaml"):
            vals = Props.read_from(file)
        else:
            err = f"File {file} not recognized as .pkl, .json, .yml or .yaml"
            raise ValueError(err)
        fk = ChannelProcKey.get_filekey_from_pattern(Path(file).name)
        out_dicts[fk.timestamp] = vals
    return out_dicts


def split_files_by_run(files):
    # sort files in dictionary where keys are first timestamp from run
    if isinstance(files, list):
        files = []
        for file in files:
            with Path(file).open() as f:
                files += f.read().splitlines()
    else:
        with Path(files).open() as f:
            files = f.read().splitlines()

    files = sorted(
        np.unique(files)
    )  # need this as sometimes files get double counted as it somehow puts in the p%-* filelist and individual runs also

    final_dict = {}
    all_file = run_splitter(sorted(files))
    for filelist in all_file:
        fk = ProcessingFileKey.get_filekey_from_pattern(Path(sorted(filelist)[0]).name)
        final_dict[fk.timestamp] = sorted(filelist)

    return final_dict, files


def save_dict_to_files(files, out_dict):
    for file in files:
        fk = ChannelProcKey.get_filekey_from_pattern(Path(file).name)
        if Path(file).suffix == ".pkl":
            with Path(file).open("wb") as w:
                pkl.dump(out_dict[fk.timestamp], w, protocol=pkl.HIGHEST_PROTOCOL)
        elif Path(file).suffix in (".json", ".yml", ".yaml"):
            Props.write_to(file, out_dict[fk.timestamp])
        else:
            err = f"File {file} not recognized as .pkl, .json, .yml or .yaml"
            raise ValueError(err)
