from __future__ import annotations

import argparse
import logging
import os
import pathlib
import pickle as pkl
import warnings

import numpy as np
import pandas as pd
from legendmeta import LegendMetadata
from legendmeta.catalog import Props
from pygama.math.distributions import gaussian
from pygama.pargen.data_cleaning import get_tcm_pulser_ids
from pygama.pargen.lq_cal import *  # noqa: F403
from pygama.pargen.lq_cal import LQCal
from pygama.pargen.utils import load_data

log = logging.getLogger(__name__)
warnings.filterwarnings(action="ignore", category=RuntimeWarning)


def get_results_dict(lq_class):
    return {
        "cal_energy_param": lq_class.cal_energy_param,
        "DEP_means": lq_class.timecorr_df.to_dict("index"),
        "rt_correction": lq_class.dt_fit_pars,
        "cut_fit_pars": lq_class.cut_fit_pars,
        "cut_value": lq_class.cut_val,
        "sfs": lq_class.low_side_sf.to_dict("index"),
    }


def fill_plot_dict(lq_class, data, plot_options, plot_dict=None):
    if plot_dict is not None:
        for key, item in plot_options.items():
            if item["options"] is not None:
                plot_dict[key] = item["function"](lq_class, data, **item["options"])
            else:
                plot_dict[key] = item["function"](lq_class, data)
    else:
        plot_dict = {}
    return plot_dict


def lq_calibration(
    data: pd.DataFrame,
    cal_dicts: dict,
    energy_param: str,
    cal_energy_param: str,
    dt_param: str,
    eres_func: callable,
    cdf: callable = gaussian,
    selection_string: str = "",
    plot_options: dict | None = None,
    debug_mode: bool = False,
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
    dt_param: string
        The drift time parameter of choice
    eres_func: callable
        The energy resolution functions
    cdf: callable
        The CDF used for the binned fitting of LQ distributions
    selection_string: string
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
    lq: LQCal class
        The LQCal object used for the LQ calibration
    """

    lq = LQCal(
        cal_dicts,
        cal_energy_param,
        dt_param,
        eres_func,
        cdf,
        selection_string,
        debug_mode=debug_mode | args.debug,
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
    return cal_dicts, get_results_dict(lq), fill_plot_dict(lq, data, plot_options), lq


argparser = argparse.ArgumentParser()
argparser.add_argument("files", help="files", nargs="*", type=str)
argparser.add_argument("--pulser_file", help="pulser_file", type=str, required=False)
argparser.add_argument("--tcm_filelist", help="tcm_filelist", type=str, required=False)

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
argparser.add_argument("--lq_results", help="lq_results", type=str)

argparser.add_argument("-d", "--debug", help="debug_mode", action="store_true")
args = argparser.parse_args()

logging.basicConfig(level=logging.DEBUG, filename=args.log, filemode="w")
logging.getLogger("numba").setLevel(logging.INFO)
logging.getLogger("parse").setLevel(logging.INFO)
logging.getLogger("lgdo").setLevel(logging.INFO)
logging.getLogger("h5py").setLevel(logging.INFO)
logging.getLogger("matplotlib").setLevel(logging.INFO)

configs = LegendMetadata(path=args.configs)
channel_dict = configs.on(args.timestamp, system=args.datatype)["snakemake_rules"][
    "pars_hit_lqcal"
]["inputs"]["lqcal_config"][args.channel]

kwarg_dict = Props.read_from(channel_dict)

ecal_dict = Props.read_from(args.ecal_file)
cal_dict = ecal_dict["pars"]["operations"]
eres_dict = ecal_dict["results"]["ecal"]

with open(args.eres_file, "rb") as o:
    object_dict = pkl.load(o)

if kwarg_dict["run_lq"] is True:
    kwarg_dict.pop("run_lq")

    cdf = eval(kwarg_dict.pop("cdf")) if "cdf" in kwarg_dict else gaussian

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
        "lq80",
        "dt_eff",
        kwarg_dict["energy_param"],
        kwarg_dict["cal_energy_param"],
        kwarg_dict["cut_field"],
    ]

    # load data in
    data, threshold_mask = load_data(
        files,
        f"{args.channel}/dsp",
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
        with open(args.tcm_filelist) as f:
            tcm_files = f.read().splitlines()
        tcm_files = sorted(np.unique(tcm_files))
        ids, mask = get_tcm_pulser_ids(
            tcm_files, args.channel, kwarg_dict.pop("pulser_multiplicity_threshold")
        )
    else:
        msg = "No pulser file or tcm filelist provided"
        raise ValueError(msg)

    data["is_pulser"] = mask[threshold_mask]

    cal_dict, out_dict, plot_dict, obj = lq_calibration(
        data,
        selection_string=f"{kwarg_dict.pop('cut_field')}&(~is_pulser)",
        cal_dicts=cal_dict,
        eres_func=eres_func,
        cdf=cdf,
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


results_dict = dict(**eres_dict, lq=out_dict)
pathlib.Path(os.path.dirname(args.hit_pars)).mkdir(parents=True, exist_ok=True)
final_hit_dict = {
    "pars": {"operations": cal_dict},
    "results": results_dict,
}
Props.write_to(args.hit_pars, final_hit_dict)

pathlib.Path(os.path.dirname(args.lq_results)).mkdir(parents=True, exist_ok=True)
final_object_dict = dict(
    **object_dict,
    lq=obj,
)
Props.write_to(args.lq_results, final_object_dict)
with open(args.lq_results, "wb") as w:
    pkl.dump(final_object_dict, w, protocol=pkl.HIGHEST_PROTOCOL)
