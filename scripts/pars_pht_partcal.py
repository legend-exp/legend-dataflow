from __future__ import annotations

import argparse
import copy
import pickle as pkl
import re
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import pygama.math.distributions as pgf
import pygama.math.histogram as pgh
from legendmeta import LegendMetadata, TextDB
from legendmeta.catalog import Props
from pygama.math.distributions import nb_poly
from pygama.pargen.data_cleaning import get_tcm_pulser_ids
from pygama.pargen.energy_cal import FWHMLinear, FWHMQuadratic, HPGeCalibration
from pygama.pargen.utils import load_data
from util.FileKey import ChannelProcKey, ProcessingFileKey
from utils.log import build_log

warnings.filterwarnings(action="ignore", category=RuntimeWarning)
warnings.filterwarnings(action="ignore", category=np.RankWarning)


def run_splitter(files):
    """
    Returns list containing lists of each run
    """

    runs = []
    run_files = []
    for file in files:
        fk = ProcessingFileKey.get_filekey_from_pattern(Path(file).name)
        if f"{fk.period}-{fk.run}" not in runs:
            runs.append(f"{fk.period}-{fk.run}")
            run_files.append([])
        for i, run in enumerate(runs):
            if run == f"{fk.period}-{fk.run}":
                run_files[i].append(file)
    return run_files


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


def bin_spectrum(
    data,
    cal_energy_param,
    selection_string,
    cut_field="is_valid_cal",
    pulser_field="is_pulser",
    erange=(0, 3000),
    dx=0.25,
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


def get_results_dict(ecal_class, data, cal_energy_param, selection_string):
    if np.isnan(ecal_class.pars).all():
        return {}
    else:
        results_dict = copy.deepcopy(ecal_class.results["hpge_fit_energy_peaks"])

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

        out_dict = {
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
            "peak_param": results_dict["peak_param"],
        }
        if "calibration_parameters" in results_dict:
            out_dict["calibration_parameters"] = results_dict["calibration_parameters"].to_dict()
            out_dict["calibration_uncertainty"] = results_dict[
                "calibration_uncertainties"
            ].to_dict()

        return out_dict


def calibrate_partition(
    data,
    cal_dicts,
    results_dicts,
    object_dicts,
    plot_dicts,
    timestamp,
    chmap,
    configs,
    channel,
    datatype,
    gen_plots=True,
):

    det_status = chmap[channel]["analysis"]["usability"]

    configs = LegendMetadata(path=configs)
    channel_dict = configs.on(timestamp, system=datatype)["snakemake_rules"]["pars_pht_partcal"][
        "inputs"
    ]["pars_pht_partcal_config"][channel]

    kwarg_dict = Props.read_from(channel_dict)

    if "plot_options" in kwarg_dict:
        for field, item in kwarg_dict["plot_options"].items():
            kwarg_dict["plot_options"][field]["function"] = eval(item["function"])

    nruns = len(np.unique(data["run_timestamp"]))

    # calibrate
    pk_pars = [
        # (238.632, (10, 10), pgf.gauss_on_step), #double line, Pb-212
        # (240.986, (10, 10), pgf.gauss_on_step), #double line, Ra-224
        (277.371, (10, 7), pgf.gauss_on_linear),  # Tl-208
        (288.2, (7, 10), pgf.gauss_on_linear),  # Bi-212
        (300.087, (10, 10), pgf.gauss_on_linear),  # Pb-212
        (452.98, (10, 10), pgf.gauss_on_linear),  # Bi-212
        # (511, (20, 20), pgf.gauss_on_step), double line, #e+e-
        (549.73, (10, 10), pgf.gauss_on_linear),  # Rn-220
        (583.187, (20, 20), pgf.hpge_peak),  # Tl-208
        (727.330, (20, 20), pgf.hpge_peak),  # Bi-212
        (763.13, (20, 10), pgf.gauss_on_linear),  # Tl-208
        (785.37, (10, 20), pgf.gauss_on_linear),  # Bi-212
        (860.557, (20, 20), pgf.hpge_peak),  # Tl-208
        (893.408, (20, 20), pgf.gauss_on_linear),  # Bi-212
        (927.6, (20, 20), pgf.gauss_on_linear),  # Tl-208
        (952.120, (20, 20), pgf.gauss_on_linear),  # Bi-212
        (982.7, (20, 20), pgf.gauss_on_linear),  # Tl-208
        (1078.62, (20, 7), pgf.gauss_on_linear),  # Bi-212
        (1093.9, (7, 20), pgf.gauss_on_linear),  # Tl-208
        (1512.7, (20, 20), pgf.gauss_on_linear),  # Bi-212
        (1592.511, (20, 20), pgf.hpge_peak),  # Tl-208 DEP
        (1620.50, (20, 20), pgf.hpge_peak),  # Bi-212
        (1679.7, (20, 20), pgf.gauss_on_linear),  # Bi-212
        (1806.0, (20, 20), pgf.gauss_on_linear),  # Bi-212
        (2103.511, (20, 20), pgf.hpge_peak),  # Tl-208 SEP
        (2614.511, (40, 20), pgf.hpge_peak),  # Tl-208
        (3125.511, (20, 20), pgf.gauss_on_linear),  # Summation
        (3197.7, (20, 20), pgf.gauss_on_linear),  # Summation
        (3475.1, (20, 20), pgf.gauss_on_linear),
    ]

    glines = [pk_par[0] for pk_par in pk_pars]

    if "cal_energy_params" not in kwarg_dict:
        cal_energy_params = [energy_param + "_cal" for energy_param in kwarg_dict["energy_params"]]
    else:
        cal_energy_params = kwarg_dict["cal_energy_params"]

    selection_string = f"~is_pulser&{kwarg_dict['final_cut_field']}"

    ecal_results = {}
    partcal_plot_dict = {}
    full_object_dict = {}

    for energy_param, cal_energy_param in zip(kwarg_dict["energy_params"], cal_energy_params):
        energy = data.query(selection_string)[energy_param].to_numpy()
        full_object_dict[cal_energy_param] = HPGeCalibration(
            energy_param,
            glines,
            1,
            kwarg_dict.get("deg", 0),
            debug_mode=kwarg_dict.get("debug_mode", False) | args.debug,  # , fixed={1: 1}
        )
        full_object_dict[cal_energy_param].hpge_get_energy_peaks(
            energy,
            etol_kev=5 if det_status == "on" else 10,
            update_cal_pars=bool(det_status == "on"),
        )
        if det_status != "on":
            full_object_dict[cal_energy_param].peak_locs = np.array(glines)

        full_object_dict[cal_energy_param].hpge_fit_energy_peaks(
            energy,
            peak_pars=pk_pars,
            tail_weight=kwarg_dict.get("tail_weight", 0),
            n_events=kwarg_dict.get("n_events", None),
            allowed_p_val=kwarg_dict.get("p_val", 0),
            update_cal_pars=bool(det_status == "on"),
            bin_width_kev=0.1 if nruns > 3 else 0.5,
        )

        if (
            2614.511 not in full_object_dict[cal_energy_param].peaks_kev
            and det_status == "on"
            and (cal_energy_param == "cuspEmax_ctc_cal")
        ):
            csqr = full_object_dict[cal_energy_param].results["hpge_fit_energy_peaks"][
                "peak_parameters"
            ][2614.511]["chi_square"]
            if csqr[0] / csqr[1] < 100:
                allowed_p_val = (
                    0.9
                    * full_object_dict[cal_energy_param].results["hpge_fit_energy_peaks"][
                        "peak_parameters"
                    ][2614.511]["p_value"]
                )

                full_object_dict[cal_energy_param] = HPGeCalibration(
                    energy_param,
                    [*full_object_dict[cal_energy_param].peaks_kev, 2614.511],
                    1,
                    kwarg_dict.get("deg", 0),  # , fixed={1: 1}
                )
                full_object_dict[cal_energy_param].hpge_get_energy_peaks(
                    energy,
                    etol_kev=5 if det_status == "on" else 10,
                    update_cal_pars=bool(det_status == "on"),
                )

                full_object_dict[cal_energy_param].hpge_fit_energy_peaks(
                    energy,
                    peak_pars=pk_pars,
                    tail_weight=kwarg_dict.get("tail_weight", 0),
                    n_events=kwarg_dict.get("n_events", None),
                    allowed_p_val=allowed_p_val,
                    update_cal_pars=bool(det_status == "on"),
                    bin_width_kev=0.2 if nruns > 3 else 0.5,
                )
            else:
                err = f"2614.511 peak not found in {cal_energy_param} fit, reduced csqr {csqr[0]/csqr[1]} not below 10, check fit"
                raise ValueError(err)

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

        ecal_results[cal_energy_param] = get_results_dict(
            full_object_dict[cal_energy_param], data, cal_energy_param, selection_string
        )
        cal_dicts = update_cal_dicts(
            cal_dicts, {cal_energy_param: full_object_dict[cal_energy_param].gen_pars_dict()}
        )
        if "ctc" in cal_energy_param:
            no_ctc_dict = full_object_dict[cal_energy_param].gen_pars_dict()
            no_ctc_dict["expression"] = no_ctc_dict["expression"].replace("ctc", "noctc")

            cal_dicts = update_cal_dicts(
                cal_dicts, {cal_energy_param.replace("ctc", "noctc"): no_ctc_dict}
            )
            cal_dicts = update_cal_dicts(
                cal_dicts,
                {
                    cal_energy_param.replace("_ctc", ""): {
                        "expression": f"where({cal_energy_param.replace('ctc', 'noctc')}>{kwarg_dict.get('dt_theshold_kev',100)}, {cal_energy_param}, {cal_energy_param.replace('ctc','noctc')})",
                        "parameters": {},
                    }
                },
            )

        if gen_plots is True:
            param_plot_dict = {}
            if ~np.isnan(full_object_dict[cal_energy_param].pars).all():
                param_plot_dict["fwhm_fit"] = full_object_dict[cal_energy_param].plot_eres_fit(
                    energy
                )
                param_plot_dict["cal_fit"] = full_object_dict[cal_energy_param].plot_cal_fit(
                    energy
                )
                if det_status == "on":
                    param_plot_dict["cal_fit_with_errors"] = full_object_dict[
                        cal_energy_param
                    ].plot_cal_fit_with_errors(energy)
                if (
                    len(
                        full_object_dict[cal_energy_param].results["hpge_fit_energy_peaks"][
                            "peak_parameters"
                        ]
                    )
                    < 17
                ):
                    param_plot_dict["peak_fits"] = full_object_dict[cal_energy_param].plot_fits(
                        energy, ncols=4, nrows=4
                    )
                elif (
                    len(
                        full_object_dict[cal_energy_param].results["hpge_fit_energy_peaks"][
                            "peak_parameters"
                        ]
                    )
                    < 26
                ):
                    param_plot_dict["peak_fits"] = full_object_dict[cal_energy_param].plot_fits(
                        energy, ncols=5, nrows=5
                    )
                else:
                    param_plot_dict["peak_fits"] = full_object_dict[cal_energy_param].plot_fits(
                        energy, ncols=6, nrows=5
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
            partcal_plot_dict[cal_energy_param] = param_plot_dict

        for peak_dict in (
            full_object_dict[cal_energy_param]
            .results["hpge_fit_energy_peaks"]["peak_parameters"]
            .values()
        ):
            peak_dict["function"] = peak_dict["function"].name
            peak_dict["parameters"] = peak_dict["parameters"].to_dict()
            peak_dict["uncertainties"] = peak_dict["uncertainties"].to_dict()

    out_result_dicts = {}
    for tstamp, result_dict in results_dicts.items():
        out_result_dicts[tstamp] = dict(**result_dict, partition_ecal=ecal_results)

    out_object_dicts = {}
    for tstamp, object_dict in object_dicts.items():
        out_object_dicts[tstamp] = dict(**object_dict, partition_ecal=full_object_dict)

    common_dict = partcal_plot_dict.pop("common") if "common" in list(partcal_plot_dict) else None
    out_plot_dicts = {}
    for tstamp, plot_dict in plot_dicts.items():
        if "common" in list(plot_dict) and common_dict is not None:
            plot_dict["common"].update(partcal_plot_dict["common"])
        elif common_dict is not None:
            plot_dict["common"] = common_dict
        plot_dict.update({"partition_ecal": partcal_plot_dict})
        out_plot_dicts[tstamp] = plot_dict

    return cal_dicts, out_result_dicts, out_object_dicts, out_plot_dicts


if __name__ == "__main__":
    argparser = argparse.ArgumentParser()
    argparser.add_argument("--input_files", help="files", type=str, nargs="*", required=True)
    argparser.add_argument(
        "--pulser_files", help="pulser_file", nargs="*", type=str, required=False
    )
    argparser.add_argument(
        "--tcm_filelist", help="tcm_filelist", type=str, nargs="*", required=False
    )
    argparser.add_argument("--ecal_file", help="ecal_file", type=str, nargs="*", required=True)
    argparser.add_argument("--eres_file", help="eres_file", type=str, nargs="*", required=True)
    argparser.add_argument("--inplots", help="eres_file", type=str, nargs="*", required=True)

    argparser.add_argument("--timestamp", help="Datatype", type=str, required=True)
    argparser.add_argument("--datatype", help="Datatype", type=str, required=True)
    argparser.add_argument("--channel", help="Channel", type=str, required=True)

    argparser.add_argument("--configs", help="configs", type=str, required=True)
    argparser.add_argument("--metadata", help="metadata path", type=str, required=True)
    argparser.add_argument("--log", help="log_file", type=str)

    argparser.add_argument("--plot_file", help="plot_file", type=str, nargs="*", required=False)
    argparser.add_argument("--hit_pars", help="hit_pars", nargs="*", type=str)
    argparser.add_argument("--fit_results", help="fit_results", nargs="*", type=str)

    argparser.add_argument("-d", "--debug", help="debug_mode", action="store_true")
    args = argparser.parse_args()

    configs = TextDB(args.configs, lazy=True).on(args.timestamp, system=args.datatype)
    config_dict = configs["snakemake_rules"]["pars_pht_partcal"]

    log = build_log(config_dict, args.log)

    meta = LegendMetadata(path=args.metadata)
    chmap = meta.channelmap(args.timestamp, system=args.datatype)
    channel = f"ch{chmap[args.channel].daq.rawid:07}"

    cal_dict = {}
    results_dicts = {}
    for ecal in args.ecal_file:
        cal = Props.read_from(ecal)

        fk = ChannelProcKey.get_filekey_from_pattern(Path(ecal).name)
        cal_dict[fk.timestamp] = cal["pars"]
        results_dicts[fk.timestamp] = cal["results"]

    object_dict = {}
    for ecal in args.eres_file:
        with Path(ecal).open("rb") as o:
            cal = pkl.load(o)
        fk = ChannelProcKey.get_filekey_from_pattern(Path(ecal).name)
        object_dict[fk.timestamp] = cal

    inplots_dict = {}
    if args.inplots:
        for ecal in args.inplots:
            with Path(ecal).open("rb") as o:
                cal = pkl.load(o)
            fk = ChannelProcKey.get_filekey_from_pattern(Path(ecal).name)
            inplots_dict[fk.timestamp] = cal

    # sort files in dictionary where keys are first timestamp from run
    files = []
    for file in args.input_files:
        with Path(file).open() as f:
            files += f.read().splitlines()

    files = sorted(
        np.unique(files)
    )  # need this as sometimes files get double counted as it somehow puts in the p%-* filelist and individual runs also

    final_dict = {}
    all_file = run_splitter(sorted(files))
    for filelist in all_file:
        fk = ProcessingFileKey.get_filekey_from_pattern(Path(sorted(filelist)[0]).name)
        timestamp = fk.timestamp
        final_dict[timestamp] = sorted(filelist)

    channel_dict = config_dict["inputs"]["pars_pht_partcal_config"][args.channel]
    kwarg_dict = Props.read_from(channel_dict)

    params = [
        kwarg_dict["final_cut_field"],
        "timestamp",
    ]
    params += kwarg_dict["energy_params"]

    # load data in
    data, threshold_mask = load_data(
        final_dict,
        f"{channel}/dsp",
        cal_dict,
        params=params,
        threshold=kwarg_dict["threshold"],
        return_selection_mask=True,
        cal_energy_param=kwarg_dict["energy_params"][0],
    )

    if args.pulser_files:
        mask = np.array([], dtype=bool)
        for file in args.pulser_files:
            pulser_dict = Props.read_from(file)
            pulser_mask = np.array(pulser_dict["mask"])
            mask = np.append(mask, pulser_mask)
        if "pulser_multiplicity_threshold" in kwarg_dict:
            kwarg_dict.pop("pulser_multiplicity_threshold")

    elif args.tcm_filelist:
        # get pulser mask from tcm files
        with Path(args.tcm_filelist).open() as f:
            tcm_files = f.read().splitlines()
        tcm_files = sorted(np.unique(tcm_files))
        ids, mask = get_tcm_pulser_ids(
            tcm_files, channel, kwarg_dict["pulser_multiplicity_threshold"]
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

    cal_dicts, results_dicts, object_dicts, plot_dicts = calibrate_partition(
        data,
        cal_dict,
        results_dicts,
        object_dict,
        inplots_dict,
        timestamp,
        chmap,
        args.configs,
        args.channel,
        args.datatype,
        gen_plots=bool(args.plot_file),
    )

    if args.plot_file:
        for plot_file in args.plot_file:
            Path(plot_file).parent.mkdir(parents=True, exist_ok=True)
            with Path(plot_file).open("wb") as w:
                pkl.dump(plot_dicts[fk.timestamp], w, protocol=pkl.HIGHEST_PROTOCOL)

    for out in sorted(args.hit_pars):
        fk = ChannelProcKey.get_filekey_from_pattern(Path(out).name)
        final_hit_dict = {
            "pars": cal_dict[fk.timestamp],
            "results": results_dicts[fk.timestamp],
        }
        Path(out).parent.mkdir(parents=True, exist_ok=True)
        Props.write_to(out, final_hit_dict)

    for out in args.fit_results:
        fk = ChannelProcKey.get_filekey_from_pattern(Path(out).name)
        Path(out).parent.mkdir(parents=True, exist_ok=True)
        with Path(out).open("wb") as w:
            pkl.dump(object_dict[fk.timestamp], w, protocol=pkl.HIGHEST_PROTOCOL)
