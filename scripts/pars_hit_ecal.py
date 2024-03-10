from __future__ import annotations

import argparse
import json
import logging
import os
import pathlib
import pickle as pkl
import warnings
from datetime import datetime

import lgdo.lh5 as lh5
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pygama.math.histogram as pgh
from legendmeta import LegendMetadata
from legendmeta.catalog import Props
from matplotlib.colors import LogNorm
from pygama.pargen.ecal_th import *  # noqa: F403
from pygama.pargen.ecal_th import apply_cuts, calibrate_parameter
from pygama.pargen.utils import get_tcm_pulser_ids, load_data
from scipy.stats import binned_statistic

log = logging.getLogger(__name__)
mpl.use("agg")
sto = lh5.LH5Store()

warnings.filterwarnings(action="ignore", category=RuntimeWarning)


def plot_baseline_timemap(
    data,
    figsize=(12, 8),
    fontsize=12,
    parameter="bl_mean",
    dx=1,
    n_spread=5,
    time_dx=180,
):
    plt.rcParams["figure.figsize"] = figsize
    plt.rcParams["font.size"] = fontsize

    time_bins = np.arange(
        (np.amin(data["timestamp"]) // time_dx) * time_dx,
        ((np.amax(data["timestamp"]) // time_dx) + 2) * time_dx,
        time_dx,
    )

    mean = np.nanpercentile(data[parameter], 50)
    spread = mean - np.nanpercentile(data[parameter], 10)
    fig = plt.figure()
    plt.hist2d(
        data["timestamp"],
        data[parameter],
        bins=[
            time_bins,
            np.arange(mean - n_spread * spread, mean + n_spread * spread + dx, dx),
        ],
        norm=LogNorm(),
    )

    ticks, labels = plt.xticks()
    plt.xlabel(f"Time starting : {datetime.utcfromtimestamp(ticks[0]).strftime('%d/%m/%y %H:%M')}")
    plt.ylabel("Baseline Value")
    plt.ylim([mean - n_spread * spread, mean + n_spread * spread])

    plt.xticks(
        ticks,
        [datetime.utcfromtimestamp(tick).strftime("%H:%M") for tick in ticks],
    )
    plt.close()
    return fig


def bin_bl_stability(data, time_slice=180, parameter="bl_mean"):
    utime_array = data["timestamp"]
    select_bls = data[parameter].to_numpy()

    time_bins = np.arange(
        (np.amin(utime_array) // time_slice) * time_slice,
        ((np.amax(utime_array) // time_slice) + 2) * time_slice,
        time_slice,
    )
    # bin time values
    times_average = (time_bins[:-1] + time_bins[1:]) / 2

    def nanmedian(x):
        return np.nanpercentile(x, 50) if len(x) >= 10 else np.nan

    def error(x):
        return np.nanvar(x) / np.sqrt(len(x)) if len(x) >= 10 else np.nan

    par_average, _, _ = binned_statistic(
        data["timestamp"], select_bls, statistic=nanmedian, bins=time_bins
    )
    par_error, _, _ = binned_statistic(
        data["timestamp"], select_bls, statistic=error, bins=time_bins
    )

    return {"time": times_average, "baseline": par_average, "spread": par_error}


def bin_baseline(data, parameter="bl_mean-baseline", dx=1, bl_range=None):
    if bl_range is None:
        bl_range = [-500, 500]
    par_array = data.eval(parameter)
    bins = np.arange(bl_range[0], bl_range[1], dx)
    bl_array, bins, _ = pgh.get_hist(par_array, bins=bins)
    return {"bl_array": bl_array, "bins": (bins[1:] + bins[:-1]) / 2}


def baseline_tracking_plots(files, lh5_path, plot_options=None):
    if plot_options is None:
        plot_options = {}
    plot_dict = {}
    data = sto.read(lh5_path, files, field_mask=["bl_mean", "baseline", "timestamp"])[0].view_as(
        "pd"
    )
    for key, item in plot_options.items():
        if item["options"] is not None:
            plot_dict[key] = item["function"](data, **item["options"])
        else:
            plot_dict[key] = item["function"](data)
    return plot_dict

def get_results_dict(ecal_class, data):
    if np.isnan(ecal_class.pars).all():
        return {}
    else:
        fwhm_linear = ecal_class.fwhm_fit_linear.copy()
        fwhm_linear["parameters"] = fwhm_linear["parameters"].to_dict()
        fwhm_linear["uncertainties"] = fwhm_linear["uncertainties"].to_dict()
        fwhm_linear["cov"] = fwhm_linear["cov"].tolist()
        fwhm_quad = ecal_class.fwhm_fit_quadratic.copy()
        fwhm_quad["parameters"] = fwhm_quad["parameters"].to_dict()
        fwhm_quad["uncertainties"] = fwhm_quad["uncertainties"].to_dict()
        fwhm_quad["cov"] = fwhm_quad["cov"].tolist()

        pk_dict = {
            Ei: {
                "function": func_i.__name__,
                "module": func_i.__module__,
                "parameters_in_ADC": parsi.to_dict(),
                "uncertainties_in_ADC": errorsi.to_dict(),
                "p_val": pvali,
                "fwhm_in_keV": list(fwhmi),
                "pk_position":(posi, posuni),
            }
            for i, (Ei, parsi, errorsi, pvali, fwhmi, posi, posuni, func_i) in enumerate(
                zip(
                    ecal_class.results["fitted_keV"],
                    ecal_class.results["pk_pars"][ecal_class.results["pk_validities"]],
                    ecal_class.results["pk_errors"][ecal_class.results["pk_validities"]],
                    ecal_class.results["pk_pvals"][ecal_class.results["pk_validities"]],
                    ecal_class.results["pk_fwhms"],
                    ecal_class.results["pk_pos"],
                    ecal_class.results["pk_pos_uncertainties"],
                    ecal_class.funcs,
                )
            )
        }

        return {
            "total_fep": len(
                data.query(
                    f"{ecal_class.cal_energy_param}>2604&{ecal_class.cal_energy_param}<2624"
                )
            ),
            "total_dep": len(
                data.query(
                    f"{ecal_class.cal_energy_param}>1587&{ecal_class.cal_energy_param}<1597"
                )
            ),
            "pass_fep": len(
                data.query(
                    f"{ecal_class.cal_energy_param}>2604&{ecal_class.cal_energy_param}<2624&{ecal_class.selection_string}"
                )
            ),
            "pass_dep": len(
                data.query(
                    f"{ecal_class.cal_energy_param}>1587&{ecal_class.cal_energy_param}<1597&{ecal_class.selection_string}"
                )
            ),
            "eres_linear": fwhm_linear,
            "eres_quadratic": fwhm_quad,
            "fitted_peaks": ecal_class.results["fitted_keV"].tolist(),
            "pk_fits": pk_dict,
            "mode":ecal_class.results["mode"],
        }

def energy_cal_th(
    data: pd.Dataframe,
    energy_params: list[str],
    cal_energy_params: list | None = None,
    selection_string: str = "",
    hit_dict: dict | None = None,
    cut_parameters: dict[str, int] | None = None,
    plot_options: dict | None = None,
    threshold: int = 0,
    p_val: float = 0,
    n_events: int | None = None,
    final_cut_field: str = "is_valid_cal",
    simplex: bool = True,
    guess_keV: float | None = None,
    tail_weight=100,
    deg: int = 1,
) -> tuple(dict, dict, dict, dict):
    data, hit_dict = apply_cuts(
        data,
        hit_dict if hit_dict is not None else {},
        cut_parameters if cut_parameters is not None else {},
        final_cut_field,
    )

    if cal_energy_params is None:
        cal_energy_params = [energy_param + "_cal" for energy_param in energy_params]


    glines = [
        # 238.632,
        583.191,
        727.330,
        860.564,
        1592.53,
        1620.50,
        2103.53,
        2614.50,
    ]  # gamma lines used for calibration
    range_keV = [
        # (8, 8),
        (20, 20),
        (30, 30),
        (30, 30),
        (40, 20),
        (20, 40),
        (40, 40),
        (60, 60),
    ]  # side bands width
    funcs = [
        # pgf.extended_gauss_step_pdf,
        pgf.extended_radford_pdf,
        pgf.extended_radford_pdf,
        pgf.extended_radford_pdf,
        pgf.extended_radford_pdf,
        pgf.extended_radford_pdf,
        pgf.extended_radford_pdf,
        pgf.extended_radford_pdf,
    ]
    gof_funcs = [
        # pgf.gauss_step_pdf,
        pgf.radford_pdf,
        pgf.radford_pdf,
        pgf.radford_pdf,
        pgf.radford_pdf,
        pgf.radford_pdf,
        pgf.radford_pdf,
        pgf.radford_pdf,
    ]

    results_dict = {}
    plot_dict = {}
    full_object_dict = {}
    for energy_param, cal_energy_param in zip(energy_params, cal_energy_params):
        full_object_dict[cal_energy_param] = calibrate_parameter(
            energy_param,
            glines,
            range_keV,
            funcs,
            gof_funcs,
            selection_string,
            plot_options,
            guess_keV,
            threshold,
            p_val,
            n_events,
            simplex,
            deg,
            tail_weight=tail_weight,
            cal_energy_param=cal_energy_param,
        )
        full_object_dict[cal_energy_param].calibrate_parameter(data)
        results_dict[cal_energy_param] = get_results_dict(full_object_dict[cal_energy_param], data)
        hit_dict.update(full_object_dict[cal_energy_param].hit_dict)
        if ~np.isnan(full_object_dict[cal_energy_param].pars).all():
            plot_dict[cal_energy_param] = (
                full_object_dict[cal_energy_param].fill_plot_dict(data).copy()
            )

    log.info("Finished all calibrations")
    return hit_dict, results_dict, plot_dict, full_object_dict


if __name__ == "__main__":
    argparser = argparse.ArgumentParser()
    argparser.add_argument("--files", help="files", nargs="*", type=str)
    argparser.add_argument("--tcm_filelist", help="tcm_filelist", type=str, required=True)
    argparser.add_argument("--ctc_dict", help="ctc_dict", nargs="*")
    argparser.add_argument("--in_hit_dict", help="in_hit_dict", nargs="*", required=False)
    argparser.add_argument("--inplot_dict", help="inplot_dict", nargs="*", required=False)

    argparser.add_argument("--configs", help="config", type=str, required=True)
    argparser.add_argument("--datatype", help="Datatype", type=str, required=True)
    argparser.add_argument("--timestamp", help="Timestamp", type=str, required=True)
    argparser.add_argument("--channel", help="Channel", type=str, required=True)
    argparser.add_argument("--tier", help="tier", type=str, default="hit")

    argparser.add_argument("--log", help="log_file", type=str)

    argparser.add_argument("--plot_path", help="plot_path", type=str, required=False)
    argparser.add_argument("--save_path", help="save_path", type=str)
    argparser.add_argument("--results_path", help="results_path", type=str)
    args = argparser.parse_args()

    logging.basicConfig(level=logging.DEBUG, filename=args.log, filemode="w")
    logging.getLogger("numba").setLevel(logging.INFO)
    logging.getLogger("parse").setLevel(logging.INFO)
    logging.getLogger("lgdo").setLevel(logging.INFO)
    logging.getLogger("h5py").setLevel(logging.INFO)
    logging.getLogger("matplotlib").setLevel(logging.INFO)
    logging.getLogger("legendmeta").setLevel(logging.INFO)

    if args.in_hit_dict:
        hit_dict = Props.read_from(args.in_hit_dict)

    database_dic = Props.read_from(args.ctc_dict)

    hit_dict = hit_dict.update(database_dic[args.channel]["ctc_params"])

    # get metadata dictionary
    configs = LegendMetadata(path=args.configs)
    channel_dict = configs.on(args.timestamp, system=args.datatype)["snakemake_rules"]
    if args.tier == "hit":
        channel_dict = channel_dict["pars_hit_ecal"]["inputs"]["ecal_config"][args.channel]
    elif args.tier == "pht":
        channel_dict = channel_dict["pars_pht_ecal"]["inputs"]["ecal_config"][args.channel]
    else:
        msg = "invalid tier"
        raise ValueError(msg)

    kwarg_dict = Props.read_from(channel_dict)

    # convert plot functions from strings to functions and split off baseline and common plots
    for field, item in kwarg_dict["plot_options"].items():
        kwarg_dict["plot_options"][field]["function"] = eval(item["function"])

    bl_plots = kwarg_dict.pop("bl_plot_options")
    for field, item in bl_plots.items():
        bl_plots[field]["function"] = eval(item["function"])
    common_plots = kwarg_dict.pop("common_plots")

    # load data in
    data, threshold_mask = load_data(
        args.files,
        f"{args.channel}/dsp",
        hit_dict,
        params=kwarg_dict["energy_params"]
        + list(kwarg_dict["cut_parameters"])
        + ["timestamp", "trapTmax"],
        threshold=kwarg_dict["threshold"],
        return_selection_mask=True,
        cal_energy_param="trapTmax",
    )

    # get pulser mask from tcm files
    with open(args.tcm_filelist) as f:
        tcm_files = f.read().splitlines()
    tcm_files = sorted(np.unique(tcm_files))
    ids, mask = get_tcm_pulser_ids(
        tcm_files, args.channel, kwarg_dict.pop("pulser_multiplicity_threshold")
    )
    data["is_pulser"] = mask[threshold_mask]

    # run energy calibration
    out_dict, result_dict, plot_dict, ecal_object = energy_cal_th(
        data,
        hit_dict=hit_dict,
        selection_string=f"({kwarg_dict['final_cut_field']})&(~is_pulser)",
        **kwarg_dict,
    )

    # get baseline plots and save all plots to file
    if args.plot_path:
        common_dict = baseline_tracking_plots(
            sorted(args.files), f"{args.channel}/dsp", plot_options=bl_plots
        )

        for plot in list(common_dict):
            if plot not in common_plots:
                plot_item = common_dict.pop(plot)
                plot_dict.update({plot: plot_item})

        for key, item in plot_dict.items():
            if isinstance(item, dict) and len(item) > 0:
                param_dict = {}
                for plot in common_plots:
                    if plot in item:
                        param_dict.update({plot: item[plot]})
                common_dict.update({key: param_dict})
        plot_dict = {"ecal":plot_dict}
        plot_dict["common"] = common_dict

        if args.inplot_dict:
            with open(args.inplot_dict, "rb") as f:
                total_plot_dict = pkl.load(args.inplot_dict, protocol=pkl.HIGHEST_PROTOCOL)
            if "common" in total_plot_dict:
                total_plot_dict["common"].update(common_dict)
            else:
                plot_dict["common"] = common_dict

            total_plot_dict = total_plot_dict.update(plot_dict)

        pathlib.Path(os.path.dirname(args.plot_path)).mkdir(parents=True, exist_ok=True)
        with open(args.plot_path, "wb") as f:
            pkl.dump(plot_dict, f, protocol=pkl.HIGHEST_PROTOCOL)

    # save output dictionary
    output_dict = {"pars": out_dict, "results": {"ecal":result_dict}}
    with open(args.save_path, "w") as fp:
        pathlib.Path(os.path.dirname(args.save_path)).mkdir(parents=True, exist_ok=True)
        json.dump(output_dict, fp, indent=4)

    # save calibration objects
    with open(args.results_path, "wb") as fp:
        pathlib.Path(os.path.dirname(args.results_path)).mkdir(parents=True, exist_ok=True)
        pkl.dump({"ecal":ecal_object}, fp, protocol=pkl.HIGHEST_PROTOCOL)
