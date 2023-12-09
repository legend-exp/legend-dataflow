from __future__ import annotations

import argparse
import json
import logging
import os
import pathlib
import pickle as pkl

import numpy as np
import pandas as pd
from legendmeta import LegendMetadata
from legendmeta.catalog import Props
from pygama.math.peak_fitting import gauss_cdf
from pygama.pargen.lq_cal import *  # noqa: F403
from pygama.pargen.lq_cal import cal_lq
from pygama.pargen.utils import get_tcm_pulser_ids, load_data
from util.FileKey import ChannelProcKey, ProcessingFileKey

log = logging.getLogger(__name__)


def lq_calibration(
    data: pd.DataFrame,
    cal_dicts: dict,
    energy_param: str,
    cal_energy_param: str,
    eres_func: callable,
    cdf: callable = gauss_cdf,
    selection_string: str = "",
    plot_options: dict | None = None,
):
    """Loads in data from the provided files and runs the LQ calibration on said files

    Parameters
    ----------
    data: pd.DataFrame
        A dataframe containing the data used for calibrating LQ
    cal_dicts: dict
        A dict of hit-level operations to apply to the data
    energy_param: string
        The energy parameter of choice. Used for normalizing the
        raw lq values
    cal_energy_param: string
        The calibrated energy parameter of choice
    eres_func: callable
        The energy resolution functions
    cdf: callable
        The CDF used for the binned fitting of LQ distributions
    cut_field: string
        A string of flags to apply to the data when running the calibration
    plot_options: dict
        A dict containing the plot functions the user wants to run,and any
        user options to provide those plot functions
    Returns
    -------
    cal_dicts: dict
        The user provided dict, updated with hit-level operations for LQ
    results_dict: dict
        A dict containing the results of the LQ calibration
    plot_dict: dict
        A dict containing all the figures specified by the plot options
    lq: cal_lq class
        The cal_lq object used for the LQ calibration
    """

    lq = cal_lq(
        cal_dicts,
        cal_energy_param,
        eres_func,
        cdf,
        selection_string,
        plot_options,
    )

    data["LQ_Ecorr"] = np.divide(data["lq80"], data[energy_param])

    lq.update_cal_dicts(
        {
            "LQ_Ecorr": {
                "expression": f"lq80/{energy_param}",
                "parameters": {},
            }
        }
    )

    lq.calibrate(data, "LQ_Ecorr")
    log.info("Calibrated LQ")
    return cal_dicts, lq.get_results_dict(), lq.fill_plot_dict(data), lq


argparser = argparse.ArgumentParser()
argparser.add_argument("--input_files", help="files", type=str, nargs="*", required=True)
argparser.add_argument("--tcm_filelist", help="tcm_filelist", type=str, nargs="*", required=True)
argparser.add_argument("--ecal_file", help="ecal_file", type=str, nargs="*", required=True)
argparser.add_argument("--eres_file", help="eres_file", type=str, nargs="*", required=True)
argparser.add_argument("--inplots", help="eres_file", type=str, nargs="*", required=True)

argparser.add_argument("--configs", help="configs", type=str, required=True)
argparser.add_argument("--timestamp", help="Datatype", type=str, required=True)
argparser.add_argument("--datatype", help="Datatype", type=str, required=True)
argparser.add_argument("--channel", help="Channel", type=str, required=True)

argparser.add_argument("--log", help="log_file", type=str)

argparser.add_argument("--plot_file", help="plot_file", type=str, nargs="*", required=False)
argparser.add_argument("--hit_pars", help="hit_pars", nargs="*", type=str)
argparser.add_argument("--lq_results", help="lq_results", nargs="*", type=str)
args = argparser.parse_args()

logging.basicConfig(level=logging.DEBUG, filename=args.log, filemode="w")
logging.getLogger("numba").setLevel(logging.INFO)
logging.getLogger("parse").setLevel(logging.INFO)
logging.getLogger("lgdo").setLevel(logging.INFO)
logging.getLogger("h5py").setLevel(logging.INFO)
logging.getLogger("matplotlib").setLevel(logging.INFO)


def run_splitter(files):
    """
    Returns list containing lists of each run
    """

    runs = []
    run_files = []
    for file in files:
        fk = ProcessingFileKey.get_filekey_from_pattern(os.path.basename(file))
        if f"{fk.period}-{fk.run}" not in runs:
            runs.append(f"{fk.period}-{fk.run}")
            run_files.append([])
        for i, run in enumerate(runs):
            if run == f"{fk.period}-{fk.run}":
                run_files[i].append(file)
    return run_files


configs = LegendMetadata(path=args.configs)
channel_dict = configs.on(args.timestamp, system=args.datatype)["snakemake_rules"][
    "pars_pht_lqcal"
]["inputs"]["lqcal_config"][args.channel]

kwarg_dict = Props.read_from(channel_dict)

cal_dict = {}
results_dicts = {}
if isinstance(args.ecal_file, list):
    for ecal in args.ecal_file:
        with open(ecal) as o:
            cal = json.load(o)

        fk = ChannelProcKey.get_filekey_from_pattern(os.path.basename(ecal))
        cal_dict[fk.timestamp] = cal["pars"]["operations"]
        results_dicts[fk.timestamp] = cal["results"]
else:
    with open(args.ecal_file) as o:
        cal = json.load(o)

    fk = ChannelProcKey.get_filekey_from_pattern(os.path.basename(args.ecal_file))
    cal_dict[fk.timestamp] = cal["pars"]["operations"]
    results_dicts[fk.timestamp] = cal["results"]

object_dict = {}
if isinstance(args.eres_file, list):
    for ecal in args.eres_file:
        with open(ecal, "rb") as o:
            cal = pkl.load(o)
        fk = ChannelProcKey.get_filekey_from_pattern(os.path.basename(ecal))
        object_dict[fk.timestamp] = cal
else:
    with open(args.eres_file, "rb") as o:
        cal = pkl.load(o)
    fk = ChannelProcKey.get_filekey_from_pattern(os.path.basename(args.eres_file))
    object_dict[fk.timestamp] = cal

inplots_dict = {}
if args.inplots:
    if isinstance(args.inplots, list):
        for ecal in args.inplots:
            with open(ecal, "rb") as o:
                cal = pkl.load(o)
            fk = ChannelProcKey.get_filekey_from_pattern(os.path.basename(ecal))
            inplots_dict[fk.timestamp] = cal
    else:
        with open(args.inplots, "rb") as o:
            cal = pkl.load(o)
        fk = ChannelProcKey.get_filekey_from_pattern(os.path.basename(args.inplots))
        inplots_dict[fk.timestamp] = cal

if "plot_options" in kwarg_dict:
    for field, item in kwarg_dict["plot_options"].items():
        kwarg_dict["plot_options"][field]["function"] = eval(item["function"])


# sort files in dictionary where keys are first timestamp from run
if isinstance(args.input_files, list):
    files = []
    for file in args.input_files:
        with open(file) as f:
            files += f.read().splitlines()
else:
    with open(args.input_files) as f:
        files = f.read().splitlines()

files = sorted(
    np.unique(files)
)  # need this as sometimes files get double counted as it somehow puts in the p%-* filelist and individual runs also

final_dict = {}
all_file = run_splitter(sorted(files))
for filelist in all_file:
    fk = ProcessingFileKey.get_filekey_from_pattern(os.path.basename(sorted(filelist)[0]))
    timestamp = fk.timestamp
    final_dict[timestamp] = sorted(filelist)

# run lq cal
if kwarg_dict.pop("run_lq") is True:
    params = [
        "lq80",
        "dt_eff",
        kwarg_dict["energy_param"],
        kwarg_dict["cal_energy_param"],
        kwarg_dict["cut_field"],
    ]

    # load data in
    data, threshold_mask = load_data(
        final_dict,
        f"{args.channel}/dsp",
        cal_dict,
        params=params,
        threshold=kwarg_dict.pop("threshold"),
        return_selection_mask=True,
    )

    # get pulser mask from tcm files
    if isinstance(args.tcm_filelist, list):
        tcm_files = []
        for file in args.tcm_filelist:
            with open(file) as f:
                tcm_files += f.read().splitlines()
    else:
        with open(args.tcm_filelist) as f:
            tcm_files = f.read().splitlines()

    tcm_files = sorted(np.unique(tcm_files))
    ids, mask = get_tcm_pulser_ids(
        tcm_files, args.channel, kwarg_dict.pop("pulser_multiplicity_threshold")
    )
    data["is_pulser"] = mask[threshold_mask]

    for tstamp in cal_dict:
        if tstamp not in np.unique(data["run_timestamp"]):
            row = {key: [False] if data.dtypes[key] == "bool" else [np.nan] for key in data}
            row["run_timestamp"] = tstamp
            row = pd.DataFrame(row)
            data = pd.concat([data, row])

    cdf = eval(kwarg_dict.pop("cdf")) if "cdf" in kwarg_dict else gauss_cdf

    try:
        eres = results_dicts[list(results_dicts)[0]]["partition_ecal"][
            kwarg_dict["cal_energy_param"]
        ]["eres_linear"].copy()

        def eres_func(x):
            return eval(eres["expression"], dict(x=x, **eres["parameters"]))

        if np.isnan(eres_func(2000)):
            raise RuntimeError
    except (KeyError, RuntimeError):
        try:
            eres = results_dicts[list(results_dicts)[0]]["ecal"][kwarg_dict["cal_energy_param"]][
                "eres_linear"
            ].copy()

            def eres_func(x):
                return eval(eres["expression"], dict(x=x, **eres["parameters"]))

        except KeyError:

            def eres_func(x):
                return x * np.nan

    cal_dict, out_dict, plot_dict, lq_obj = lq_calibration(
        data,
        selection_string=f"{kwarg_dict.pop('cut_field')}&(~is_pulser)",
        cal_dicts=cal_dict,
        eres_func=eres_func,
        cdf=cdf,
        **kwarg_dict,
    )

    # need to change eres func as can't pickle lambdas
    try:
        lq_obj.eres_func = results_dicts[list(results_dicts)[0]]["partition_ecal"][
            kwarg_dict["cal_energy_param"]
        ]["eres_linear"].copy()
    except KeyError:
        lq_obj.eres_func = {}
else:
    out_dict = {}
    plot_dict = {}
    lq_obj = None


if args.plot_file:
    common_dict = plot_dict.pop("common") if "common" in list(plot_dict) else None

    if isinstance(args.plot_file, list):
        for plot_file in args.plot_file:
            fk = ChannelProcKey.get_filekey_from_pattern(os.path.basename(plot_file))
            if args.inplots:
                out_plot_dict = inplots_dict[fk.timestamp]
                out_plot_dict.update({"lq": plot_dict})
            else:
                out_plot_dict = {"lq": plot_dict}
            if "common" in list(out_plot_dict) and common_dict is not None:
                out_plot_dict["common"].update(common_dict)
            elif common_dict is not None:
                out_plot_dict["common"] = common_dict
            pathlib.Path(os.path.dirname(plot_file)).mkdir(parents=True, exist_ok=True)
            with open(plot_file, "wb") as w:
                pkl.dump(out_plot_dict, w, protocol=pkl.HIGHEST_PROTOCOL)
    else:
        if args.inplots:
            fk = ChannelProcKey.get_filekey_from_pattern(os.path.basename(args.plot_file))
            out_plot_dict = inplots_dict[fk.timestamp]
            out_plot_dict.update({"lq": plot_dict})
        else:
            out_plot_dict = {"lq": plot_dict}
        if "common" in list(out_plot_dict) and common_dict is not None:
            out_plot_dict["common"].update(common_dict)
        elif common_dict is not None:
            out_plot_dict["common"] = common_dict
        pathlib.Path(os.path.dirname(args.plot_file)).mkdir(parents=True, exist_ok=True)
        with open(args.plot_file, "wb") as w:
            pkl.dump(out_plot_dict, w, protocol=pkl.HIGHEST_PROTOCOL)


for out in sorted(args.hit_pars):
    fk = ChannelProcKey.get_filekey_from_pattern(os.path.basename(out))
    final_hit_dict = {
        "pars": {"operations": cal_dict[fk.timestamp]},
        "results": dict(
            **results_dicts[fk.timestamp],
            lq=out_dict,
        ),
    }
    pathlib.Path(os.path.dirname(out)).mkdir(parents=True, exist_ok=True)
    with open(out, "w") as w:
        json.dump(final_hit_dict, w, indent=4)

for out in args.lq_results:
    fk = ChannelProcKey.get_filekey_from_pattern(os.path.basename(out))
    final_object_dict = dict(
        **object_dict[fk.timestamp],
        lq=lq_obj,
    )
    pathlib.Path(os.path.dirname(out)).mkdir(parents=True, exist_ok=True)
    with open(out, "wb") as w:
        pkl.dump(final_object_dict, w, protocol=pkl.HIGHEST_PROTOCOL)
