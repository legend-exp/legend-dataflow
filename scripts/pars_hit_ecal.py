from __future__ import annotations

import argparse
import copy
import json
import logging
import os
import pathlib
import pickle as pkl
import warnings
from datetime import datetime

os.environ["PYGAMA_PARALLEL"] = "false"
os.environ["PYGAMA_FASTMATH"] = "false"

import lgdo.lh5 as lh5
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pygama.math.distributions as pgf
import pygama.math.histogram as pgh
from legendmeta import LegendMetadata
from legendmeta.catalog import Props
from matplotlib.colors import LogNorm
from pygama.math.distributions import nb_poly
from pygama.pargen.data_cleaning import get_mode_stdev, get_tcm_pulser_ids
from pygama.pargen.energy_cal import FWHMLinear, FWHMQuadratic, HPGeCalibration
from pygama.pargen.utils import load_data
from scipy.stats import binned_statistic

log = logging.getLogger(__name__)
mpl.use("agg")
sto = lh5.LH5Store()

warnings.filterwarnings(action="ignore", category=RuntimeWarning)
warnings.filterwarnings(action="ignore", category=np.RankWarning)


def plot_2614_timemap(
    data,
    cal_energy_param,
    selection_string,
    figsize=(12, 8),
    fontsize=12,
    erange=(2580, 2630),
    dx=1,
    time_dx=180,
):
    plt.rcParams["figure.figsize"] = figsize
    plt.rcParams["font.size"] = fontsize

    selection = data.query(f"{cal_energy_param}>2560&{cal_energy_param}<2660&{selection_string}")

    fig = plt.figure()
    if len(selection) == 0:
        pass
    else:
        time_bins = np.arange(
            (np.amin(data["timestamp"]) // time_dx) * time_dx,
            ((np.amax(data["timestamp"]) // time_dx) + 2) * time_dx,
            time_dx,
        )

        plt.hist2d(
            selection["timestamp"],
            selection[cal_energy_param],
            bins=[time_bins, np.arange(erange[0], erange[1] + dx, dx)],
            norm=LogNorm(),
        )

    ticks, labels = plt.xticks()
    plt.xlabel(f"Time starting : {datetime.utcfromtimestamp(ticks[0]).strftime('%d/%m/%y %H:%M')}")
    plt.ylabel("Energy(keV)")
    plt.ylim([erange[0], erange[1]])

    plt.xticks(
        ticks,
        [datetime.utcfromtimestamp(tick).strftime("%H:%M") for tick in ticks],
    )
    plt.close()
    return fig


def plot_pulser_timemap(
    data,
    cal_energy_param,
    selection_string,  # noqa: ARG001
    pulser_field="is_pulser",
    figsize=(12, 8),
    fontsize=12,
    dx=0.2,
    time_dx=180,
    n_spread=3,
):
    plt.rcParams["figure.figsize"] = figsize
    plt.rcParams["font.size"] = fontsize

    time_bins = np.arange(
        (np.amin(data["timestamp"]) // time_dx) * time_dx,
        ((np.amax(data["timestamp"]) // time_dx) + 2) * time_dx,
        time_dx,
    )

    selection = data.query(pulser_field)
    fig = plt.figure()
    if len(selection) == 0:
        pass

    else:
        mean = np.nanpercentile(selection[cal_energy_param], 50)
        spread = mean - np.nanpercentile(selection[cal_energy_param], 10)

        plt.hist2d(
            selection["timestamp"],
            selection[cal_energy_param],
            bins=[
                time_bins,
                np.arange(mean - n_spread * spread, mean + n_spread * spread + dx, dx),
            ],
            norm=LogNorm(),
        )
        plt.ylim([mean - n_spread * spread, mean + n_spread * spread])
    ticks, labels = plt.xticks()
    plt.xlabel(f"Time starting : {datetime.utcfromtimestamp(ticks[0]).strftime('%d/%m/%y %H:%M')}")
    plt.ylabel("Energy(keV)")

    plt.xticks(
        ticks,
        [datetime.utcfromtimestamp(tick).strftime("%H:%M") for tick in ticks],
    )
    plt.close()
    return fig


def get_median(x):
    if len(x[~np.isnan(x)]) >= 10:
        return np.nan
    else:
        return np.nanpercentile(x, 50)


def get_err(x):
    if len(x[~np.isnan(x)]) >= 10:
        return np.nan
    else:
        return np.nanvar(x) / np.sqrt(len(x))


def bin_pulser_stability(
    data,
    cal_energy_param,
    selection_string,  # noqa: ARG001
    pulser_field="is_pulser",
    time_slice=180,
):
    selection = data.query(pulser_field)

    utime_array = data["timestamp"]
    select_energies = selection[cal_energy_param].to_numpy()

    time_bins = np.arange(
        (np.amin(utime_array) // time_slice) * time_slice,
        ((np.amax(utime_array) // time_slice) + 2) * time_slice,
        time_slice,
    )
    # bin time values
    times_average = (time_bins[:-1] + time_bins[1:]) / 2

    if len(selection) == 0:
        return {
            "time": times_average,
            "energy": np.full_like(times_average, np.nan),
            "spread": np.full_like(times_average, np.nan),
        }

    par_average, _, _ = binned_statistic(
        selection["timestamp"], select_energies, statistic=get_median, bins=time_bins
    )
    par_error, _, _ = binned_statistic(
        selection["timestamp"], select_energies, statistic=get_err, bins=time_bins
    )

    return {"time": times_average, "energy": par_average, "spread": par_error}


def bin_stability(
    data,
    cal_energy_param,
    selection_string,
    time_slice=180,
    energy_range=(2585, 2660),
):
    selection = data.query(
        f"{cal_energy_param}>{energy_range[0]}&{cal_energy_param}<{energy_range[1]}&{selection_string}"
    )

    utime_array = data["timestamp"]
    select_energies = selection[cal_energy_param].to_numpy()

    time_bins = np.arange(
        (np.amin(utime_array) // time_slice) * time_slice,
        ((np.amax(utime_array) // time_slice) + 2) * time_slice,
        time_slice,
    )
    # bin time values
    times_average = (time_bins[:-1] + time_bins[1:]) / 2

    if len(selection) == 0:
        return {
            "time": times_average,
            "energy": np.full_like(times_average, np.nan),
            "spread": np.full_like(times_average, np.nan),
        }

    par_average, _, _ = binned_statistic(
        selection["timestamp"], select_energies, statistic=get_median, bins=time_bins
    )
    par_error, _, _ = binned_statistic(
        selection["timestamp"], select_energies, statistic=get_err, bins=time_bins
    )

    return {"time": times_average, "energy": par_average, "spread": par_error}


def bin_spectrum(
    data,
    cal_energy_param,
    selection_string,
    cut_field="is_valid_cal",
    pulser_field="is_pulser",
    erange=(0, 3000),
    dx=2,
):
    bins = np.arange(erange[0], erange[1] + dx, dx)
    return {
        "bins": pgh.get_bin_centers(bins),
        "counts": np.histogram(data.query(selection_string)[cal_energy_param], bins)[0],
        "cut_counts": np.histogram(
            data.query(f"(~{cut_field})&(~{pulser_field})")[cal_energy_param],
            bins,
        )[0],
        "pulser_counts": np.histogram(
            data.query(pulser_field)[cal_energy_param],
            bins,
        )[0],
    }


def bin_survival_fraction(
    data,
    cal_energy_param,
    selection_string,
    cut_field="is_valid_cal",
    pulser_field="is_pulser",
    erange=(0, 3000),
    dx=6,
):
    counts_pass, bins_pass, _ = pgh.get_hist(
        data.query(selection_string)[cal_energy_param],
        bins=np.arange(erange[0], erange[1] + dx, dx),
    )
    counts_fail, bins_fail, _ = pgh.get_hist(
        data.query(f"(~{cut_field})&(~{pulser_field})")[cal_energy_param],
        bins=np.arange(erange[0], erange[1] + dx, dx),
    )
    sf = 100 * (counts_pass + 10 ** (-6)) / (counts_pass + counts_fail + 10 ** (-6))
    return {"bins": pgh.get_bin_centers(bins_pass), "sf": sf}


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
    data = lh5.read_as(lh5_path, files, "pd", field_mask=["bl_mean", "baseline", "timestamp"])
    for key, item in plot_options.items():
        if item["options"] is not None:
            plot_dict[key] = item["function"](data, **item["options"])
        else:
            plot_dict[key] = item["function"](data)
    return plot_dict


def monitor_parameters(files, lh5_path, parameters):
    data = lh5.read_as(lh5_path, files, "pd", field_mask=parameters)
    out_dict = {}
    for param in parameters:
        mode, stdev = get_mode_stdev(data[param].to_numpy())
        out_dict[param] = {"mode": mode, "stdev": stdev}
    return out_dict


def get_results_dict(ecal_class, data, cal_energy_param, selection_string):
    if np.isnan(ecal_class.pars).all():
        return {}
    else:
        results_dict = copy.deepcopy(ecal_class.results["hpge_fit_energy_peaks_1"])

        if "FWHMLinear" in results_dict:
            fwhm_linear = results_dict["FWHMLinear"]
            fwhm_linear["function"] = fwhm_linear["function"].__name__
            fwhm_linear["parameters"] = fwhm_linear["parameters"].to_dict()
            fwhm_linear["uncertainties"] = fwhm_linear["uncertainties"].to_dict()
            fwhm_linear["cov"] = fwhm_linear["cov"].tolist()
        else:
            fwhm_linear = None

        if "FWHMQuadratic" in results_dict:
            fwhm_quad = results_dict["FWHMQuadratic"]
            fwhm_quad["function"] = fwhm_quad["function"].__name__
            fwhm_quad["parameters"] = fwhm_quad["parameters"].to_dict()
            fwhm_quad["uncertainties"] = fwhm_quad["uncertainties"].to_dict()
            fwhm_quad["cov"] = fwhm_quad["cov"].tolist()
        else:
            fwhm_quad = None

        pk_dict = results_dict["peak_parameters"]

        for _, dic in pk_dict.items():
            dic["function"] = dic["function"].name
            dic["parameters"] = dic["parameters"].to_dict()
            dic["uncertainties"] = dic["uncertainties"].to_dict()
            dic.pop("covariance")

        return {
            "total_fep": len(data.query(f"{cal_energy_param}>2604&{cal_energy_param}<2624")),
            "total_dep": len(data.query(f"{cal_energy_param}>1587&{cal_energy_param}<1597")),
            "pass_fep": len(
                data.query(f"{cal_energy_param}>2604&{cal_energy_param}<2624&{selection_string}")
            ),
            "pass_dep": len(
                data.query(f"{cal_energy_param}>1587&{cal_energy_param}<1597&{selection_string}")
            ),
            "eres_linear": fwhm_linear,
            "eres_quadratic": fwhm_quad,
            "fitted_peaks": ecal_class.peaks_kev.tolist(),
            "pk_fits": pk_dict,
        }


if __name__ == "__main__":
    argparser = argparse.ArgumentParser()
    argparser.add_argument("--files", help="filelist", nargs="*", type=str)
    argparser.add_argument("--tcm_filelist", help="tcm_filelist", type=str, required=False)
    argparser.add_argument("--pulser_file", help="pulser_file", type=str, required=False)

    argparser.add_argument("--ctc_dict", help="ctc_dict", nargs="*")
    argparser.add_argument("--in_hit_dict", help="in_hit_dict", required=False)
    argparser.add_argument("--inplot_dict", help="inplot_dict", required=False)

    argparser.add_argument("--configs", help="config", type=str, required=True)
    argparser.add_argument("--datatype", help="Datatype", type=str, required=True)
    argparser.add_argument("--timestamp", help="Timestamp", type=str, required=True)
    argparser.add_argument("--channel", help="Channel", type=str, required=True)
    argparser.add_argument("--tier", help="tier", type=str, default="hit")

    argparser.add_argument("--metadata", help="metadata path", type=str, required=True)

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

    meta = LegendMetadata(path=args.metadata)
    chmap = meta.channelmap(args.timestamp)

    det_status = chmap.map("daq.rawid")[int(args.channel[2:])]["analysis"]["usability"]

    if args.in_hit_dict:
        hit_dict = Props.read_from(args.in_hit_dict)

    database_dic = Props.read_from(args.ctc_dict)

    hit_dict.update(database_dic[args.channel]["ctc_params"])

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

    with open(args.files[0]) as f:
        files = f.read().splitlines()
    files = sorted(files)

    # load data in
    data, threshold_mask = load_data(
        files,
        f"{args.channel}/dsp",
        hit_dict,
        params=[*kwarg_dict["energy_params"], kwarg_dict["cut_param"], "timestamp", "trapTmax"],
        threshold=kwarg_dict["threshold"],
        return_selection_mask=True,
        cal_energy_param="trapTmax",
    )

    if args.pulser_file:
        with open(args.pulser_file) as f:
            pulser_dict = json.load(f)
        mask = np.array(pulser_dict["mask"])

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

    pk_pars = [
        (583.191, (20, 20), pgf.hpge_peak),
        (727.330, (30, 30), pgf.hpge_peak),
        (860.564, (30, 25), pgf.hpge_peak),
        (1592.53, (40, 20), pgf.gauss_on_step),
        (1620.50, (20, 40), pgf.gauss_on_step),
        (2103.53, (40, 40), pgf.gauss_on_step),
        (2614.50, (60, 60), pgf.hpge_peak),
    ]

    glines = [pk_par[0] for pk_par in pk_pars]

    if "cal_energy_params" not in kwarg_dict:
        cal_energy_params = [energy_param + "_cal" for energy_param in kwarg_dict["energy_params"]]
    else:
        cal_energy_params = kwarg_dict["cal_energy_params"]

    selection_string = f"~is_pulser&{kwarg_dict['cut_param']}"

    results_dict = {}
    plot_dict = {}
    full_object_dict = {}

    for energy_param, cal_energy_param in zip(kwarg_dict["energy_params"], cal_energy_params):
        e_uncal = data.query(selection_string)[energy_param].to_numpy()

        hist, bins, bar = pgh.get_hist(
            e_uncal[
                (e_uncal > np.nanpercentile(e_uncal, 95))
                & (e_uncal < np.nanpercentile(e_uncal, 99.9))
            ],
            dx=9,
            range=[np.nanpercentile(e_uncal, 95), np.nanpercentile(e_uncal, 99.9)],
        )

        guess = 2614.553 / bins[np.nanargmax(hist)]
        full_object_dict[cal_energy_param] = HPGeCalibration(
            energy_param,
            glines,
            guess,
            kwarg_dict.get("deg", 0),
        )
        full_object_dict[cal_energy_param].hpge_get_energy_peaks(
            e_uncal, etol_kev=5 if det_status == "on" else 10
        )
        got_peaks_kev = full_object_dict[cal_energy_param].peaks_kev.copy()
        full_object_dict[cal_energy_param].hpge_fit_energy_peaks(
            e_uncal,
            peaks_kev=[2614.50],
            peak_pars=pk_pars,
            tail_weight=kwarg_dict.get("tail_weight", 0),
            n_events=kwarg_dict.get("n_events", None),
            allowed_p_val=kwarg_dict.get("p_val", 0),
            update_cal_pars=bool(det_status == "on"),
        )
        full_object_dict[cal_energy_param].hpge_fit_energy_peaks(
            e_uncal,
            peaks_kev=got_peaks_kev,
            peak_pars=pk_pars,
            tail_weight=kwarg_dict.get("tail_weight", 0),
            n_events=kwarg_dict.get("n_events", None),
            allowed_p_val=kwarg_dict.get("p_val", 0),
            update_cal_pars=False,
        )

        full_object_dict[cal_energy_param].get_energy_res_curve(
            FWHMLinear,
            interp_energy_kev={"Qbb": 2039.0},
        )
        full_object_dict[cal_energy_param].get_energy_res_curve(
            FWHMQuadratic,
            interp_energy_kev={"Qbb": 2039.0},
        )

        data[cal_energy_param] = nb_poly(
            data[energy_param].to_numpy(), full_object_dict[cal_energy_param].pars
        )

        results_dict[cal_energy_param] = get_results_dict(
            full_object_dict[cal_energy_param], data, cal_energy_param, selection_string
        )

        hit_dict.update({cal_energy_param: full_object_dict[cal_energy_param].gen_pars_dict()})
        if args.plot_path:
            param_plot_dict = {}
            if ~np.isnan(full_object_dict[cal_energy_param].pars).all():
                param_plot_dict["fwhm_fit"] = full_object_dict[cal_energy_param].plot_eres_fit(
                    e_uncal
                )
                param_plot_dict["cal_fit"] = full_object_dict[cal_energy_param].plot_cal_fit(
                    e_uncal
                )
                param_plot_dict["peak_fits"] = full_object_dict[cal_energy_param].plot_fits(
                    e_uncal
                )

                if "plot_options" in kwarg_dict:
                    for key, item in kwarg_dict["plot_options"].items():
                        if item["options"] is not None:
                            param_plot_dict[key] = item["function"](
                                data,
                                cal_energy_param,
                                selection_string,
                                **item["options"],
                            )
                        else:
                            param_plot_dict[key] = item["function"](
                                data,
                                cal_energy_param,
                                selection_string,
                            )
            plot_dict[cal_energy_param] = param_plot_dict

        for peak_dict in (
            full_object_dict[cal_energy_param]
            .results["hpge_fit_energy_peaks_1"]["peak_parameters"]
            .values()
        ):
            peak_dict["function"] = peak_dict["function"].name
            peak_dict["parameters"] = peak_dict["parameters"].to_dict()
            peak_dict["uncertainties"] = peak_dict["uncertainties"].to_dict()
        for peak_dict in (
            full_object_dict[cal_energy_param]
            .results["hpge_fit_energy_peaks"]["peak_parameters"]
            .values()
        ):
            peak_dict["function"] = peak_dict["function"].name
            peak_dict["parameters"] = peak_dict["parameters"].to_dict()
            peak_dict["uncertainties"] = peak_dict["uncertainties"].to_dict()

    if "monitoring_parameters" in kwarg_dict:
        monitor_dict = monitor_parameters(
            files, f"{args.channel}/dsp", kwarg_dict["monitoring_parameters"]
        )
        results_dict.update({"monitoring_parameters": monitor_dict})

    # get baseline plots and save all plots to file
    if args.plot_path:
        common_dict = baseline_tracking_plots(
            sorted(files), f"{args.channel}/dsp", plot_options=bl_plots
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

        if args.inplot_dict:
            with open(args.inplot_dict, "rb") as f:
                total_plot_dict = pkl.load(f)
        else:
            total_plot_dict = {}

        if "common" in total_plot_dict:
            total_plot_dict["common"].update(common_dict)
        else:
            total_plot_dict["common"] = common_dict

        total_plot_dict.update({"ecal": plot_dict})

        pathlib.Path(os.path.dirname(args.plot_path)).mkdir(parents=True, exist_ok=True)
        with open(args.plot_path, "wb") as f:
            pkl.dump(total_plot_dict, f, protocol=pkl.HIGHEST_PROTOCOL)

    # save output dictionary
    output_dict = {"pars": hit_dict, "results": {"ecal": results_dict}}
    with open(args.save_path, "w") as fp:
        pathlib.Path(os.path.dirname(args.save_path)).mkdir(parents=True, exist_ok=True)
        json.dump(output_dict, fp, indent=4)

    # save calibration objects
    with open(args.results_path, "wb") as fp:
        pathlib.Path(os.path.dirname(args.results_path)).mkdir(parents=True, exist_ok=True)
        pkl.dump({"ecal": full_object_dict}, fp, protocol=pkl.HIGHEST_PROTOCOL)
