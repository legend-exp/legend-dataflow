from __future__ import annotations

import argparse
import copy
import json
import logging
import os
import pathlib
import pickle as pkl
import warnings
from typing import Callable

os.environ["PYGAMA_PARALLEL"] = "false"
os.environ["PYGAMA_FASTMATH"] = "false"

import numpy as np
import pandas as pd
from legendmeta import LegendMetadata
from legendmeta.catalog import Props
from pygama.pargen.AoE_cal import *  # noqa: F403
from pygama.pargen.AoE_cal import CalAoE, Pol1, SigmaFit, aoe_peak
from pygama.pargen.data_cleaning import get_tcm_pulser_ids
from pygama.pargen.utils import load_data
from util.FileKey import ChannelProcKey, ProcessingFileKey

log = logging.getLogger(__name__)
warnings.filterwarnings(action="ignore", category=RuntimeWarning)


def get_results_dict(aoe_class):
    result_dict = {}
    for tstamp in aoe_class.low_side_sfs_by_run:
        result_dict[tstamp] = {
            "cal_energy_param": aoe_class.cal_energy_param,
            "dt_param": aoe_class.dt_param,
            "rt_correction": aoe_class.dt_corr,
            "1000-1300keV": aoe_class.timecorr_df.to_dict("index"),
            "correction_fit_results": aoe_class.energy_corr_res_dict,
            "low_cut": aoe_class.low_cut_val,
            "high_cut": aoe_class.high_cut_val,
            "low_side_sfs": aoe_class.low_side_sfs.to_dict("index"),
            "2_side_sfs": aoe_class.two_side_sfs.to_dict("index"),
            "low_side_sfs_by_run": aoe_class.low_side_sfs_by_run[tstamp].to_dict("index"),
            "2_side_sfs_by_run": aoe_class.two_side_sfs_by_run[tstamp].to_dict("index"),
        }
    return result_dict


def fill_plot_dict(aoe_class, data, plot_options, plot_dict=None):
    if plot_dict is None:
        plot_dict = {}
    for key, item in plot_options.items():
        if item["options"] is not None:
            plot_dict[key] = item["function"](aoe_class, data, **item["options"])
        else:
            plot_dict[key] = item["function"](aoe_class, data)

    return plot_dict


def aoe_calibration(
    data: pd.Dataframe,
    cal_dicts: dict,
    current_param: str,
    energy_param: str,
    cal_energy_param: str,
    eres_func: Callable,
    pdf: Callable = aoe_peak,
    selection_string: str = "",
    dt_corr: bool = False,
    dep_correct: bool = False,
    dt_cut: dict | None = None,
    high_cut_val: int = 3,
    mean_func: Callable = Pol1,
    sigma_func: Callable = SigmaFit,
    # dep_acc: float = 0.9,
    dt_param: str = "dt_eff",
    comptBands_width: int = 20,
    plot_options: dict | None = None,
):
    data["AoE_Uncorr"] = data[current_param] / data[energy_param]
    aoe = CalAoE(
        cal_dicts=cal_dicts,
        cal_energy_param=cal_energy_param,
        eres_func=eres_func,
        pdf=pdf,
        selection_string=selection_string,
        dt_corr=dt_corr,
        dep_correct=dep_correct,
        dt_cut=dt_cut,
        dt_param=dt_param,
        high_cut_val=high_cut_val,
        mean_func=mean_func,
        sigma_func=sigma_func,
        compt_bands_width=comptBands_width,
    )
    aoe.update_cal_dicts(
        {
            "AoE_Uncorr": {
                "expression": f"{current_param}/{energy_param}",
                "parameters": {},
            }
        }
    )
    aoe.calibrate(data, "AoE_Uncorr")
    log.info("Calibrated A/E")
    return cal_dicts, get_results_dict(aoe), fill_plot_dict(aoe, data, plot_options), aoe


argparser = argparse.ArgumentParser()
argparser.add_argument("--input_files", help="files", type=str, nargs="*", required=True)
argparser.add_argument("--pulser_files", help="pulser_file", nargs="*", type=str, required=False)
argparser.add_argument("--tcm_filelist", help="tcm_filelist", type=str, nargs="*", required=False)
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
argparser.add_argument("--aoe_results", help="aoe_results", nargs="*", type=str)
args = argparser.parse_args()

logging.basicConfig(level=logging.DEBUG, filename=args.log, filemode="w")
logging.getLogger("numba").setLevel(logging.INFO)
logging.getLogger("parse").setLevel(logging.INFO)
logging.getLogger("lgdo").setLevel(logging.INFO)
logging.getLogger("h5py").setLevel(logging.INFO)
logging.getLogger("matplotlib").setLevel(logging.INFO)
logging.getLogger("legendmeta").setLevel(logging.INFO)


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
    "pars_pht_aoecal"
]["inputs"]["par_pht_aoecal_config"][args.channel]

kwarg_dict = Props.read_from(channel_dict)

cal_dict = {}
results_dicts = {}
if isinstance(args.ecal_file, list):
    for ecal in args.ecal_file:
        cal = Props.read_from(ecal)

        fk = ChannelProcKey.get_filekey_from_pattern(os.path.basename(ecal))
        cal_dict[fk.timestamp] = cal["pars"]
        results_dicts[fk.timestamp] = cal["results"]
else:
    cal = Props.read_from(args.ecal_file)

    fk = ChannelProcKey.get_filekey_from_pattern(os.path.basename(args.ecal_file))
    cal_dict[fk.timestamp] = cal["pars"]
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

# run aoe cal
if kwarg_dict.pop("run_aoe") is True:
    params = [
        kwarg_dict["final_cut_field"],
        kwarg_dict["current_param"],
        "tp_0_est",
        "tp_99",
        kwarg_dict["energy_param"],
        kwarg_dict["cal_energy_param"],
        "timestamp",
    ]
    if "dt_param" in kwarg_dict:
        params.append(kwarg_dict["dt_param"])
    else:
        params.append("dt_eff")

    if "dt_cut" in kwarg_dict and kwarg_dict["dt_cut"] is not None:
        params.append(kwarg_dict["dt_cut"]["out_param"])
        for tstamp in cal_dict:
            cal_dict[tstamp].update(kwarg_dict["dt_cut"]["cut"])

    # load data in
    data, threshold_mask = load_data(
        final_dict,
        f"{args.channel}/dsp",
        cal_dict,
        params=params,
        threshold=kwarg_dict.pop("threshold"),
        return_selection_mask=True,
    )

    if args.pulser_files:
        mask = np.array([], dtype=bool)
        for file in args.pulser_files:
            with open(file) as f:
                pulser_dict = json.load(f)
            pulser_mask = np.array(pulser_dict["mask"])
            mask = np.append(mask, pulser_mask)
        if "pulser_multiplicity_threshold" in kwarg_dict:
            kwarg_dict.pop("pulser_multiplicity_threshold")

    elif args.tcm_filelist:
        # get pulser mask from tcm files
        with open(args.tcm_filelist) as f:
            tcm_files = f.read().splitlines()
        tcm_files = sorted(np.unique(tcm_files))
        ids, mask = get_tcm_pulser_ids(
            tcm_files, args.channel, kwarg_dict["pulser_multiplicity_threshold"]
        )
    else:
        msg = "No pulser file or tcm filelist provided"
        raise ValueError(msg)

    data["is_pulser"] = mask[threshold_mask]

    for tstamp in cal_dict:
        if tstamp not in np.unique(data["run_timestamp"]):
            row = {key: [False] if data.dtypes[key] == "bool" else [np.nan] for key in data}
            row["run_timestamp"] = tstamp
            row = pd.DataFrame(row)
            data = pd.concat([data, row])

    pdf = eval(kwarg_dict.pop("pdf")) if "pdf" in kwarg_dict else aoe_peak

    mean_func = eval(kwarg_dict.pop("mean_func")) if "mean_func" in kwarg_dict else Pol1

    sigma_func = eval(kwarg_dict.pop("sigma_func")) if "sigma_func" in kwarg_dict else SigmaFit

    try:
        eres = copy.deepcopy(
            results_dicts[next(iter(results_dicts))]["partition_ecal"][
                kwarg_dict["cal_energy_param"]
            ]["eres_linear"]
        )

        def eres_func(x):
            return eval(eres["expression"], dict(x=x, **eres["parameters"]))

        if np.isnan(eres_func(2000)):
            raise RuntimeError
    except (KeyError, RuntimeError):
        try:
            eres = copy.deepcopy(
                results_dicts[next(iter(results_dicts))]["ecal"][kwarg_dict["cal_energy_param"]][
                    "eres_linear"
                ]
            )

            def eres_func(x):
                return eval(eres["expression"], dict(x=x, **eres["parameters"]))

        except KeyError:

            def eres_func(x):
                return x * np.nan

    cal_dict, out_dict, plot_dict, aoe_obj = aoe_calibration(
        data,
        selection_string=f"{kwarg_dict.pop('final_cut_field')}&(~is_pulser)",
        cal_dicts=cal_dict,
        eres_func=eres_func,
        pdf=pdf,
        mean_func=mean_func,
        sigma_func=sigma_func,
        **kwarg_dict,
    )
    aoe_obj.pdf = aoe_obj.pdf.name
    # need to change eres func as can't pickle lambdas
    try:
        aoe_obj.eres_func = results_dicts[next(iter(results_dicts))]["partition_ecal"][
            kwarg_dict["cal_energy_param"]
        ]["eres_linear"]
    except KeyError:
        aoe_obj.eres_func = {}
else:
    out_dict = {tstamp: None for tstamp in cal_dict}
    plot_dict = {}
    aoe_obj = None


if args.plot_file:
    common_dict = plot_dict.pop("common") if "common" in list(plot_dict) else None

    if isinstance(args.plot_file, list):
        for plot_file in args.plot_file:
            fk = ChannelProcKey.get_filekey_from_pattern(os.path.basename(plot_file))
            if args.inplots:
                out_plot_dict = inplots_dict[fk.timestamp]
                out_plot_dict.update({"aoe": plot_dict})
            else:
                out_plot_dict = {"aoe": plot_dict}
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
            out_plot_dict.update({"aoe": plot_dict})
        else:
            out_plot_dict = {"aoe": plot_dict}
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
            aoe=out_dict[fk.timestamp],
        ),
    }
    pathlib.Path(os.path.dirname(out)).mkdir(parents=True, exist_ok=True)
    with open(out, "w") as w:
        json.dump(final_hit_dict, w, indent=4)

for out in args.aoe_results:
    fk = ChannelProcKey.get_filekey_from_pattern(os.path.basename(out))
    final_object_dict = dict(
        **object_dict[fk.timestamp],
        aoe=aoe_obj,
    )
    pathlib.Path(os.path.dirname(out)).mkdir(parents=True, exist_ok=True)
    with open(out, "wb") as w:
        pkl.dump(final_object_dict, w, protocol=pkl.HIGHEST_PROTOCOL)
