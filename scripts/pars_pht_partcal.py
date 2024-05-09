from __future__ import annotations

import argparse
import copy
import json
import logging
import os
import pathlib
import pickle as pkl
import re
import warnings

os.environ["PYGAMA_PARALLEL"] = "false"
os.environ["PYGAMA_FASTMATH"] = "false"

import numpy as np
import pandas as pd
import pygama.math.distributions as pgf
import pygama.math.histogram as pgh
from legendmeta import LegendMetadata
from legendmeta.catalog import Props
from pygama.math.distributions import nb_poly
from pygama.pargen.data_cleaning import get_tcm_pulser_ids
from pygama.pargen.energy_cal import FWHMLinear, FWHMQuadratic, HPGeCalibration
from pygama.pargen.utils import load_data
from util.FileKey import ChannelProcKey, ProcessingFileKey

log = logging.getLogger(__name__)
warnings.filterwarnings(action="ignore", category=RuntimeWarning)


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
            "peak_param": results_dict["peak_param"],
        }


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

    argparser.add_argument("--configs", help="configs", type=str, required=True)
    argparser.add_argument("--timestamp", help="Datatype", type=str, required=True)
    argparser.add_argument("--datatype", help="Datatype", type=str, required=True)
    argparser.add_argument("--channel", help="Channel", type=str, required=True)

    argparser.add_argument("--log", help="log_file", type=str)
    argparser.add_argument("--metadata", help="metadata path", type=str, required=True)

    argparser.add_argument("--plot_file", help="plot_file", type=str, nargs="*", required=False)
    argparser.add_argument("--hit_pars", help="hit_pars", nargs="*", type=str)
    argparser.add_argument("--fit_results", help="fit_results", nargs="*", type=str)
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

    configs = LegendMetadata(path=args.configs)
    channel_dict = configs.on(args.timestamp, system=args.datatype)["snakemake_rules"][
        "pars_pht_partcal"
    ]["inputs"]["pars_pht_partcal_config"][args.channel]

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

    params = [
        kwarg_dict["final_cut_field"],
        "timestamp",
    ]
    params += kwarg_dict["energy_params"]

    # load data in
    data, threshold_mask = load_data(
        final_dict,
        f"{args.channel}/dsp",
        cal_dict,
        params=params,
        threshold=kwarg_dict["threshold"],
        return_selection_mask=True,
        cal_energy_param=kwarg_dict["energy_params"][0],
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

    pk_pars = [
        (238.632, (10, 10), pgf.gauss_on_step),
        (511, (30, 30), pgf.gauss_on_step),
        (583.191, (30, 30), pgf.hpge_peak),
        (727.330, (30, 30), pgf.hpge_peak),
        (763, (30, 15), pgf.gauss_on_step),
        (785, (15, 30), pgf.gauss_on_step),
        (860.564, (30, 25), pgf.hpge_peak),
        (893, (25, 30), pgf.gauss_on_step),
        (1079, (30, 30), pgf.gauss_on_step),
        (1513, (30, 30), pgf.gauss_on_step),
        (1592.53, (30, 20), pgf.hpge_peak),
        (1620.50, (20, 30), pgf.hpge_peak),
        (2103.53, (30, 30), pgf.hpge_peak),
        (2614.553, (30, 30), pgf.hpge_peak),
        (3125, (30, 30), pgf.gauss_on_step),
        (3198, (30, 30), pgf.gauss_on_step),
        (3474, (30, 30), pgf.gauss_on_step),
    ]

    glines = [pk_par[0] for pk_par in pk_pars]

    if "cal_energy_params" not in kwarg_dict:
        cal_energy_params = [energy_param + "_cal" for energy_param in kwarg_dict["energy_params"]]
    else:
        cal_energy_params = kwarg_dict["cal_energy_params"]

    selection_string = f"~is_pulser&{kwarg_dict['final_cut_field']}"

    ecal_results = {}
    plot_dict = {}
    full_object_dict = {}

    for energy_param, cal_energy_param in zip(kwarg_dict["energy_params"], cal_energy_params):
        energy = data.query(selection_string)[energy_param].to_numpy()
        full_object_dict[cal_energy_param] = HPGeCalibration(
            energy_param, glines, 1, kwarg_dict.get("deg", 0)  # , fixed={1: 1}
        )
        full_object_dict[cal_energy_param].hpge_get_energy_peaks(
            energy, etol_kev=5 if det_status == "on" else 10
        )

        if det_status != "on":
            full_object_dict[cal_energy_param].hpge_cal_energy_peak_tops(
                energy,
                update_cal_pars=True,
                allowed_p_val=0,
            )

        full_object_dict[cal_energy_param].hpge_fit_energy_peaks(
            energy,
            peak_pars=pk_pars,
            tail_weight=kwarg_dict.get("tail_weight", 0),
            n_events=kwarg_dict.get("n_events", None),
            allowed_p_val=kwarg_dict.get("p_val", 0),
            update_cal_pars=bool(det_status == "on"),
            bin_width_kev=0.25,
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

        ecal_results[cal_energy_param] = get_results_dict(
            full_object_dict[cal_energy_param], data, cal_energy_param, selection_string
        )
        cal_dict = update_cal_dicts(
            cal_dict, {cal_energy_param: full_object_dict[cal_energy_param].gen_pars_dict()}
        )
        if "ctc" in cal_energy_param:
            no_ctc_dict = full_object_dict[cal_energy_param].gen_pars_dict()
            no_ctc_dict["expression"] = no_ctc_dict["expression"].replace("ctc", "noctc")

            cal_dict = update_cal_dicts(
                cal_dict, {cal_energy_param.replace("ctc", "noctc"): no_ctc_dict}
            )
            cal_dict = update_cal_dicts(
                cal_dict,
                {
                    cal_energy_param.replace("_ctc", ""): {
                        "expression": f"where({cal_energy_param}>{kwarg_dict.get('dt_theshold_kev',100)}, {cal_energy_param}, {cal_energy_param.replace('ctc','noctc')})",
                        "parameters": {},
                    }
                },
            )

        if args.plot_file:
            param_plot_dict = {}
            if ~np.isnan(full_object_dict[cal_energy_param].pars).all():
                param_plot_dict["fwhm_fit"] = full_object_dict[cal_energy_param].plot_eres_fit(
                    energy
                )
                param_plot_dict["cal_fit"] = full_object_dict[cal_energy_param].plot_cal_fit(
                    energy
                )
                param_plot_dict["peak_fits"] = full_object_dict[cal_energy_param].plot_fits(
                    energy, ncols=4, nrows=5
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
            .results["hpge_fit_energy_peaks"]["peak_parameters"]
            .values()
        ):
            peak_dict["function"] = peak_dict["function"].name
            peak_dict["parameters"] = peak_dict["parameters"].to_dict()
            peak_dict["uncertainties"] = peak_dict["uncertainties"].to_dict()

        if det_status != "on":
            for peak_dict in (
                full_object_dict[cal_energy_param]
                .results["hpge_cal_energy_peak_tops"]["peak_parameters"]
                .values()
            ):
                peak_dict["function"] = peak_dict["function"].name
                peak_dict["parameters"] = peak_dict["parameters"].to_dict()
                peak_dict["uncertainties"] = peak_dict["uncertainties"].to_dict()

    if args.plot_file:
        common_dict = plot_dict.pop("common") if "common" in list(plot_dict) else None

        if isinstance(args.plot_file, list):
            for plot_file in args.plot_file:
                fk = ChannelProcKey.get_filekey_from_pattern(os.path.basename(plot_file))
                if args.inplots:
                    out_plot_dict = inplots_dict[fk.timestamp]
                    out_plot_dict.update({"partition_ecal": plot_dict})
                else:
                    out_plot_dict = {"partition_ecal": plot_dict}

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
                out_plot_dict.update({"partition_ecal": plot_dict})
            else:
                out_plot_dict = {"partition_ecal": plot_dict}
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
            "pars": cal_dict[fk.timestamp],
            "results": dict(**results_dicts[fk.timestamp], partition_ecal=ecal_results),
        }
        pathlib.Path(os.path.dirname(out)).mkdir(parents=True, exist_ok=True)
        with open(out, "w") as w:
            json.dump(final_hit_dict, w, indent=4)

    for out in args.fit_results:
        fk = ChannelProcKey.get_filekey_from_pattern(os.path.basename(out))
        final_object_dict = dict(**object_dict[fk.timestamp], partition_ecal=full_object_dict)
        pathlib.Path(os.path.dirname(out)).mkdir(parents=True, exist_ok=True)
        with open(out, "wb") as w:
            pkl.dump(final_object_dict, w, protocol=pkl.HIGHEST_PROTOCOL)
