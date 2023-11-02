from __future__ import annotations

import argparse
import json
import logging
import os
import pathlib
import pickle as pkl
from typing import Callable

import numpy as np
import pandas as pd
from legendmeta import LegendMetadata
from legendmeta.catalog import Props
from pygama.pargen.AoE_cal import *  # noqa: F403
from pygama.pargen.AoE_cal import cal_aoe, pol1, sigma_fit, standard_aoe
from pygama.pargen.utils import get_tcm_pulser_ids, load_data

log = logging.getLogger(__name__)


def aoe_calibration(
    data: pd.Dataframe,
    cal_dicts: dict,
    current_param: str,
    energy_param: str,
    cal_energy_param: str,
    eres_func: Callable,
    pdf: Callable = standard_aoe,
    selection_string: str = "",
    dt_corr: bool = False,
    dep_correct: bool = False,
    dt_cut: dict | None = None,
    high_cut_val: int = 3,
    mean_func: Callable = pol1,
    sigma_func: Callable = sigma_fit,
    dep_acc: float = 0.9,
    dt_param: str = "dt_eff",
    comptBands_width: int = 20,
    plot_options: dict | None = None,
):
    data["AoE_Uncorr"] = data[current_param] / data[energy_param]
    aoe = cal_aoe(
        cal_dicts,
        cal_energy_param,
        eres_func,
        pdf,
        selection_string,
        dt_corr,
        dep_acc,
        dep_correct,
        dt_cut,
        dt_param,
        high_cut_val,
        mean_func,
        sigma_func,
        comptBands_width,
        plot_options if plot_options is not None else {},
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
    return cal_dicts, aoe.get_results_dict(), aoe.fill_plot_dict(data), aoe


argparser = argparse.ArgumentParser()
argparser.add_argument("files", help="files", nargs="*", type=str)
argparser.add_argument("--tcm_filelist", help="tcm_filelist", type=str, required=True)
argparser.add_argument("--ecal_file", help="ecal_file", type=str, required=True)
argparser.add_argument("--eres_file", help="eres_file", type=str, required=True)
argparser.add_argument("--inplots", help="in_plot_path", type=str, required=False)

argparser.add_argument("--configs", help="configs", type=str, required=True)
argparser.add_argument("--datatype", help="Datatype", type=str, required=True)
argparser.add_argument("--timestamp", help="Timestamp", type=str, required=True)
argparser.add_argument("--channel", help="Channel", type=str, required=True)

argparser.add_argument("--log", help="log_file", type=str)

argparser.add_argument("--plot_file", help="plot_file", type=str, required=False)
argparser.add_argument("--hit_pars", help="hit_pars", type=str)
argparser.add_argument("--aoe_results", help="aoe_results", type=str)
args = argparser.parse_args()

logging.basicConfig(level=logging.DEBUG, filename=args.log, filemode="w")
logging.getLogger("numba").setLevel(logging.INFO)
logging.getLogger("parse").setLevel(logging.INFO)
logging.getLogger("lgdo").setLevel(logging.INFO)
logging.getLogger("h5py").setLevel(logging.INFO)
logging.getLogger("matplotlib").setLevel(logging.INFO)

configs = LegendMetadata(path=args.configs)
channel_dict = configs.on(args.timestamp, system=args.datatype)["snakemake_rules"][
    "pars_hit_aoecal"
]["inputs"]["aoecal_config"][args.channel]

kwarg_dict = Props.read_from(channel_dict)

with open(args.ecal_file) as o:
    ecal_dict = json.load(o)
cal_dict = ecal_dict["pars"]
eres_dict = ecal_dict["results"]

with open(args.eres_file, "rb") as o:
    object_dict = pkl.load(o)

if kwarg_dict["run_aoe"] is True:
    kwarg_dict.pop("run_aoe")

    pdf = eval(kwarg_dict.pop("pdf")) if "pdf" in kwarg_dict else standard_aoe

    sigma_func = eval(kwarg_dict.pop("sigma_func")) if "sigma_func" in kwarg_dict else sigma_fit

    mean_func = eval(kwarg_dict.pop("mean_func")) if "mean_func" in kwarg_dict else pol1

    if "plot_options" in kwarg_dict:
        for field, item in kwarg_dict["plot_options"].items():
            kwarg_dict["plot_options"][field]["function"] = eval(item["function"])

    with open(args.files[0]) as f:
        files = f.read().splitlines()
    files = sorted(files)

    try:
        eres = eres_dict[kwarg_dict["cal_energy_param"]]["eres_linear"].copy()

        def eres_func(x):
            return eval(eres["expression"], dict(x=x, **eres["parameters"]))

    except KeyError:

        def eres_func(x):
            return x * np.nan

    params = [
        kwarg_dict["current_param"],
        "tp_0_est",
        "tp_99",
        kwarg_dict["energy_param"],
        kwarg_dict["cal_energy_param"],
        kwarg_dict["cut_field"],
        "timestamp",
    ]

    if "dt_param" in kwarg_dict:
        params += kwarg_dict["dt_param"]
    else:
        params += "dt_eff"

    if "dt_cut" in kwarg_dict and kwarg_dict["dt_cut"] is not None:
        params.append(kwarg_dict["dt_cut"]["out_param"])

    # load data in
    data, threshold_mask = load_data(
        files,
        f"{args.channel}/dsp",
        cal_dict,
        params=params,
        threshold=kwarg_dict.pop("threshold"),
        return_selection_mask=True,
    )

    # get pulser mask from tcm files
    with open(args.tcm_filelist) as f:
        tcm_files = f.read().splitlines()
    tcm_files = sorted(tcm_files)
    ids, mask = get_tcm_pulser_ids(
        tcm_files, args.channel, kwarg_dict.pop("pulser_multiplicity_threshold")
    )
    data["is_pulser"] = mask[threshold_mask]

    cal_dict, out_dict, plot_dict, obj = aoe_calibration(
        data,
        cal_dicts=cal_dict,
        eres_func=eres_func,
        selection_string=f"{kwarg_dict.pop('cut_field')}&(~is_pulser)",
        pdf=pdf,
        mean_func=mean_func,
        sigma_func=sigma_func,
        **kwarg_dict,
    )

    # need to change eres func as can't pickle lambdas
    try:
        obj.eres_func = eres_dict[kwarg_dict["cal_energy_param"]]["eres_linear"].copy()
    except KeyError:
        obj.eres_func = {}
else:
    out_dict = {}
    plot_dict = {}
    obj = None

if args.plot_file:
    common_dict = plot_dict.pop("common") if "common" in list(plot_dict) else None
    if args.inplots:
        with open(args.inplots, "rb") as r:
            out_plot_dict = pkl.load(r)
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

pathlib.Path(os.path.dirname(args.hit_pars)).mkdir(parents=True, exist_ok=True)
with open(args.hit_pars, "w") as w:
    final_hit_dict = {
        "pars": {"operations": cal_dict},
        "results": {"ecal": eres_dict, "aoe": out_dict},
    }
    json.dump(final_hit_dict, w, indent=4)

pathlib.Path(os.path.dirname(args.aoe_results)).mkdir(parents=True, exist_ok=True)
with open(args.aoe_results, "wb") as w:
    pkl.dump({"ecal": object_dict, "aoe": obj}, w, protocol=pkl.HIGHEST_PROTOCOL)
