from __future__ import annotations

import argparse
import pickle as pkl
import warnings
from pathlib import Path
from typing import Callable

import numpy as np
import pandas as pd
from dbetto import TextDB
from dbetto.catalog import Props
from legenddataflow.convert_np import convert_dict_np_to_float
from legenddataflow.log import build_log
from legendmeta import LegendMetadata
from pygama.pargen.AoE_cal import *  # noqa: F403
from pygama.pargen.AoE_cal import CalAoE, Pol1, SigmaFit, aoe_peak
from pygama.pargen.data_cleaning import get_tcm_pulser_ids
from pygama.pargen.utils import load_data

warnings.filterwarnings(action="ignore", category=RuntimeWarning)


def get_results_dict(aoe_class):
    return {
        "cal_energy_param": aoe_class.cal_energy_param,
        "dt_param": aoe_class.dt_param,
        "rt_correction": aoe_class.dt_corr,
        "1000-1300keV": aoe_class.timecorr_df.to_dict("index"),
        "correction_fit_results": aoe_class.energy_corr_res_dict,
        "low_cut": aoe_class.low_cut_val,
        "high_cut": aoe_class.high_cut_val,
        "low_side_sfs": aoe_class.low_side_sfs.to_dict("index"),
        "2_side_sfs": aoe_class.two_side_sfs.to_dict("index"),
    }


def fill_plot_dict(aoe_class, data, plot_options, plot_dict=None):
    if plot_dict is not None:
        for key, item in plot_options.items():
            if item["options"] is not None:
                plot_dict[key] = item["function"](aoe_class, data, **item["options"])
            else:
                plot_dict[key] = item["function"](aoe_class, data)
    else:
        plot_dict = {}
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
    debug_mode: bool = False,
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
        debug_mode=debug_mode | args.debug,
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
argparser.add_argument("files", help="files", nargs="*", type=str)
argparser.add_argument("--pulser_file", help="pulser_file", type=str, required=False)
argparser.add_argument("--tcm_filelist", help="tcm_filelist", type=str, required=False)

argparser.add_argument("--ecal_file", help="ecal_file", type=str, required=True)
argparser.add_argument("--eres_file", help="eres_file", type=str, required=True)
argparser.add_argument("--inplots", help="in_plot_path", type=str, required=False)

argparser.add_argument("--configs", help="configs", type=str, required=True)
argparser.add_argument("--log", help="log_file", type=str)
argparser.add_argument("--metadata", help="metadata", type=str, required=True)


argparser.add_argument("--datatype", help="Datatype", type=str, required=True)
argparser.add_argument("--timestamp", help="Timestamp", type=str, required=True)
argparser.add_argument("--channel", help="Channel", type=str, required=True)

argparser.add_argument("--plot_file", help="plot_file", type=str, required=False)
argparser.add_argument("--hit_pars", help="hit_pars", type=str)
argparser.add_argument("--aoe_results", help="aoe_results", type=str)

argparser.add_argument("-d", "--debug", help="debug_mode", action="store_true")
args = argparser.parse_args()

configs = TextDB(args.configs, lazy=True).on(args.timestamp, system=args.datatype)
config_dict = configs["snakemake_rules"]["pars_hit_aoecal"]

log = build_log(config_dict, args.log)

meta = LegendMetadata(path=args.metadata)
channel_dict = meta.channelmap(args.timestamp, system=args.datatype)
channel = f"ch{channel_dict[args.channel].daq.rawid:07}"

channel_dict = config_dict["inputs"]["aoecal_config"][args.channel]
kwarg_dict = Props.read_from(channel_dict)


ecal_dict = Props.read_from(args.ecal_file)
cal_dict = ecal_dict["pars"]
eres_dict = ecal_dict["results"]["ecal"]

with Path(args.eres_file).open("rb") as o:
    object_dict = pkl.load(o)

if kwarg_dict["run_aoe"] is True:
    kwarg_dict.pop("run_aoe")

    pdf = eval(kwarg_dict.pop("pdf")) if "pdf" in kwarg_dict else aoe_peak

    sigma_func = eval(kwarg_dict.pop("sigma_func")) if "sigma_func" in kwarg_dict else SigmaFit

    mean_func = eval(kwarg_dict.pop("mean_func")) if "mean_func" in kwarg_dict else Pol1

    if "plot_options" in kwarg_dict:
        for field, item in kwarg_dict["plot_options"].items():
            kwarg_dict["plot_options"][field]["function"] = eval(item["function"])

    with Path(args.files[0]).open() as f:
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
        cal_dict.update(kwarg_dict["dt_cut"]["cut"])
        params.append(kwarg_dict["dt_cut"]["out_param"])

    # load data in
    data, threshold_mask = load_data(
        files,
        f"{channel}/dsp",
        cal_dict,
        params=params,
        threshold=kwarg_dict.pop("threshold"),
        return_selection_mask=True,
    )

    if args.pulser_file:
        pulser_dict = Props.read_from(args.pulser_file)
        mask = np.array(pulser_dict["mask"])
        if "pulser_multiplicity_threshold" in kwarg_dict:
            kwarg_dict.pop("pulser_multiplicity_threshold")

    elif args.tcm_filelist:
        # get pulser mask from tcm files
        with Path(args.tcm_filelist).open() as f:
            tcm_files = f.read().splitlines()
        tcm_files = sorted(np.unique(tcm_files))
        ids, mask = get_tcm_pulser_ids(
            tcm_files, channel, kwarg_dict.pop("pulser_multiplicity_threshold")
        )
    else:
        msg = "No pulser file or tcm filelist provided"
        raise ValueError(msg)

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
    obj.pdf = obj.pdf.name

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
        with Path(args.inplots).open("rb") as r:
            out_plot_dict = pkl.load(r)
        out_plot_dict.update({"aoe": plot_dict})
    else:
        out_plot_dict = {"aoe": plot_dict}

    if "common" in list(out_plot_dict) and common_dict is not None:
        out_plot_dict["common"].update(common_dict)
    elif common_dict is not None:
        out_plot_dict["common"] = common_dict

    Path(args.plot_file).parent.mkdir(parents=True, exist_ok=True)
    with Path(args.plot_file).open("wb") as w:
        pkl.dump(out_plot_dict, w, protocol=pkl.HIGHEST_PROTOCOL)

Path(args.hit_pars).parent.mkdir(parents=True, exist_ok=True)
results_dict = dict(**ecal_dict["results"], aoe=out_dict)
final_hit_dict = {
    "pars": {"operations": cal_dict},
    "results": results_dict,
}

final_hit_dict = convert_dict_np_to_float(final_hit_dict)

Props.write_to(args.hit_pars, final_hit_dict)

Path(args.aoe_results).parent.mkdir(parents=True, exist_ok=True)
final_object_dict = dict(
    **object_dict,
    aoe=obj,
)
with Path(args.aoe_results).open("wb") as w:
    pkl.dump(final_object_dict, w, protocol=pkl.HIGHEST_PROTOCOL)
