from __future__ import annotations

import argparse
import copy
import warnings

import numpy as np
import pandas as pd
import pygama.math.distributions as pgf
import pygama.math.histogram as pgh
from dbetto import TextDB
from dbetto.catalog import Props
from pygama.math.distributions import nb_poly
from pygama.pargen.energy_cal import FWHMLinear, FWHMQuadratic, HPGeCalibration
from pygama.pargen.utils import load_data

from .....log import build_log
from ....pulser_removal import get_pulser_mask
from .util import get_run_dict, save_dict_to_files, split_files_by_run, update_cal_dicts

warnings.filterwarnings(action="ignore", category=RuntimeWarning)
warnings.filterwarnings(action="ignore", category=np.RankWarning)


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
            "total_fep": len(
                data.query(f"{cal_energy_param}>2604&{cal_energy_param}<2624")
            ),
            "total_dep": len(
                data.query(f"{cal_energy_param}>1587&{cal_energy_param}<1597")
            ),
            "pass_fep": len(
                data.query(
                    f"{cal_energy_param}>2604&{cal_energy_param}<2624&{selection_string}"
                )
            ),
            "pass_dep": len(
                data.query(
                    f"{cal_energy_param}>1587&{cal_energy_param}<1597&{selection_string}"
                )
            ),
            "eres_linear": fwhm_linear,
            "eres_quadratic": fwhm_quad,
            "fitted_peaks": ecal_class.peaks_kev.tolist(),
            "pk_fits": pk_dict,
            "peak_param": results_dict["peak_param"],
        }
        if "calibration_parameters" in results_dict:
            out_dict["calibration_parameters"] = results_dict[
                "calibration_parameters"
            ].to_dict()
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
    channel,
    chmap,
    configs,
    gen_plots=True,
    debug_mode=False,
):
    det_status = chmap[channel]["analysis"]["usability"]

    if isinstance(configs, str | list):
        configs = Props.read_from(configs)

    if "plot_options" in configs:
        for field, item in configs["plot_options"].items():
            configs["plot_options"][field]["function"] = eval(item["function"])

    nruns = len(np.unique(data["run_timestamp"]))

    # calibrate
    pk_pars = [
        # (238.632, (10, 10), pgf.gauss_on_step), #double line, Pb-212
        # (240.986, (10, 10), pgf.gauss_on_step), #double line, Ra-224
        # enable these when understand non linearity more
        # (277.371, (10, 7), pgf.gauss_on_linear),  # Tl-208
        # (288.2, (7, 10), pgf.gauss_on_linear),  # Bi-212
        # (300.087, (10, 10), pgf.gauss_on_linear),  # Pb-212
        # (452.98, (10, 10), pgf.gauss_on_linear),  # Bi-212
        # # (511, (20, 20), pgf.gauss_on_step), double line, #e+e-
        # (549.73, (10, 10), pgf.gauss_on_linear),  # Rn-220
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

    if "cal_energy_params" not in configs:
        cal_energy_params = [
            energy_param + "_cal" for energy_param in configs["energy_params"]
        ]
    else:
        cal_energy_params = configs["cal_energy_params"]

    selection_string = f"~is_pulser&{configs['final_cut_field']}"

    ecal_results = {}
    partcal_plot_dict = {}
    full_object_dict = {}

    for energy_param, cal_energy_param in zip(
        configs["energy_params"], cal_energy_params
    ):
        energy = data.query(selection_string)[energy_param].to_numpy()
        full_object_dict[cal_energy_param] = HPGeCalibration(
            energy_param,
            glines,
            1,
            configs.get("deg", 0),
            debug_mode=configs.get("debug_mode", False) | debug_mode,  # , fixed={1: 1}
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
            tail_weight=configs.get("tail_weight", 0),
            n_events=configs.get("n_events", None),
            allowed_p_val=configs.get("p_val", 0),
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
                    * full_object_dict[cal_energy_param].results[
                        "hpge_fit_energy_peaks"
                    ]["peak_parameters"][2614.511]["p_value"]
                )

                full_object_dict[cal_energy_param] = HPGeCalibration(
                    energy_param,
                    [*full_object_dict[cal_energy_param].peaks_kev, 2614.511],
                    1,
                    configs.get("deg", 0),  # , fixed={1: 1}
                )
                full_object_dict[cal_energy_param].hpge_get_energy_peaks(
                    energy,
                    etol_kev=5 if det_status == "on" else 10,
                    update_cal_pars=bool(det_status == "on"),
                )

                full_object_dict[cal_energy_param].hpge_fit_energy_peaks(
                    energy,
                    peak_pars=pk_pars,
                    tail_weight=configs.get("tail_weight", 0),
                    n_events=configs.get("n_events", None),
                    allowed_p_val=allowed_p_val,
                    update_cal_pars=bool(det_status == "on"),
                    bin_width_kev=0.2 if nruns > 3 else 0.5,
                )
            else:
                err = f"2614.511 peak not found in {cal_energy_param} fit, reduced csqr {csqr[0] / csqr[1]} not below 10, check fit"
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
            cal_dicts,
            {cal_energy_param: full_object_dict[cal_energy_param].gen_pars_dict()},
        )
        if "ctc" in cal_energy_param:
            no_ctc_dict = full_object_dict[cal_energy_param].gen_pars_dict()
            no_ctc_dict["expression"] = no_ctc_dict["expression"].replace(
                "ctc", "noctc"
            )

            cal_dicts = update_cal_dicts(
                cal_dicts, {cal_energy_param.replace("ctc", "noctc"): no_ctc_dict}
            )
            cal_dicts = update_cal_dicts(
                cal_dicts,
                {
                    cal_energy_param.replace("_ctc", ""): {
                        "expression": f"where({cal_energy_param.replace('ctc', 'noctc')}>{configs.get('dt_theshold_kev', 100)}, {cal_energy_param}, {cal_energy_param.replace('ctc', 'noctc')})",
                        "parameters": {},
                    }
                },
            )

        if gen_plots is True:
            param_plot_dict = {}
            if ~np.isnan(full_object_dict[cal_energy_param].pars).all():
                param_plot_dict["fwhm_fit"] = full_object_dict[
                    cal_energy_param
                ].plot_eres_fit(energy)
                param_plot_dict["cal_fit"] = full_object_dict[
                    cal_energy_param
                ].plot_cal_fit(energy)
                if det_status == "on":
                    param_plot_dict["cal_fit_with_errors"] = full_object_dict[
                        cal_energy_param
                    ].plot_cal_fit_with_errors(energy)
                if (
                    len(
                        full_object_dict[cal_energy_param].results[
                            "hpge_fit_energy_peaks"
                        ]["peak_parameters"]
                    )
                    < 17
                ):
                    param_plot_dict["peak_fits"] = full_object_dict[
                        cal_energy_param
                    ].plot_fits(energy, ncols=4, nrows=4)
                elif (
                    len(
                        full_object_dict[cal_energy_param].results[
                            "hpge_fit_energy_peaks"
                        ]["peak_parameters"]
                    )
                    < 26
                ):
                    param_plot_dict["peak_fits"] = full_object_dict[
                        cal_energy_param
                    ].plot_fits(energy, ncols=5, nrows=5)
                else:
                    param_plot_dict["peak_fits"] = full_object_dict[
                        cal_energy_param
                    ].plot_fits(energy, ncols=6, nrows=5)

                if "plot_options" in configs:
                    for key, item in configs["plot_options"].items():
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

    out_result_dicts = {
        tstamp: dict(**result_dict, partition_ecal=ecal_results)
        for tstamp, result_dict in results_dicts.items()
    }
    out_object_dicts = {
        tstamp: dict(**object_dict, partition_ecal=full_object_dict)
        for tstamp, object_dict in object_dicts.items()
    }

    common_dict = (
        partcal_plot_dict.pop("common") if "common" in list(partcal_plot_dict) else None
    )
    out_plot_dicts = {}
    for tstamp, plot_dict in plot_dicts.items():
        if "common" in list(plot_dict) and common_dict is not None:
            plot_dict["common"].update(partcal_plot_dict["common"])
        elif common_dict is not None:
            plot_dict["common"] = common_dict
        plot_dict.update({"partition_ecal": partcal_plot_dict})
        out_plot_dicts[tstamp] = plot_dict

    return cal_dicts, out_result_dicts, out_object_dicts, out_plot_dicts


def par_geds_pht_ecal_part() -> None:
    argparser = argparse.ArgumentParser()
    argparser.add_argument(
        "--input-files", help="files", type=str, nargs="*", required=True
    )
    argparser.add_argument(
        "--pulser-files", help="pulser_file", nargs="*", type=str, required=False
    )
    argparser.add_argument(
        "--ecal-file", help="ecal_file", type=str, nargs="*", required=True
    )
    argparser.add_argument(
        "--eres-file", help="eres_file", type=str, nargs="*", required=True
    )
    argparser.add_argument(
        "--inplots", help="eres_file", type=str, nargs="*", required=True
    )

    argparser.add_argument("--timestamp", help="Datatype", type=str, required=True)
    argparser.add_argument("--datatype", help="Datatype", type=str, required=True)
    argparser.add_argument("--channel", help="Channel", type=str, required=True)
    argparser.add_argument("--table-name", help="table name", type=str, required=True)

    argparser.add_argument("--configs", help="configs", type=str, required=True)
    argparser.add_argument("--metadata", help="metadata path", type=str, required=True)
    argparser.add_argument("--log", help="log_file", type=str)

    argparser.add_argument(
        "--plot-file", help="plot_file", type=str, nargs="*", required=False
    )
    argparser.add_argument("--hit-pars", help="hit_pars", nargs="*", type=str)
    argparser.add_argument("--fit-results", help="fit_results", nargs="*", type=str)

    argparser.add_argument("-d", "--debug", help="debug_mode", action="store_true")
    args = argparser.parse_args()

    configs = TextDB(args.configs, lazy=True).on(args.timestamp, system=args.datatype)
    config_dict = configs["snakemake_rules"]["pars_pht_partcal"]

    build_log(config_dict, args.log)

    chmap = TextDB(path=args.metadata, lazy=True).on(
        args.timestamp, system=args.datatype
    )

    # par files
    par_dict = get_run_dict(args.ecal_file)
    cal_dicts = {tstamp: val["pars"] for tstamp, val in par_dict.items()}
    results_dicts = {tstamp: val["results"] for tstamp, val in par_dict.items()}

    # obj files
    object_dicts = get_run_dict(args.eres_file)

    # plot files
    inplots_dict = {}
    if args.inplots:
        inplots_dict = get_run_dict(args.inplots)

    # get files split by run
    final_dict, _ = split_files_by_run(args.files)

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
        args.table_name,
        cal_dicts,
        params=params,
        threshold=kwarg_dict["threshold"],
        return_selection_mask=True,
        cal_energy_param=kwarg_dict["energy_params"][0],
    )

    mask = get_pulser_mask(pulser_file=args.pulser_files)
    if "pulser_multiplicity_threshold" in kwarg_dict:
        kwarg_dict.pop("pulser_multiplicity_threshold")

    data["is_pulser"] = mask[threshold_mask]

    for tstamp in cal_dicts:
        if tstamp not in np.unique(data["run_timestamp"]):
            row = {
                key: [False] if data.dtypes[key] == "bool" else [np.nan] for key in data
            }
            row["run_timestamp"] = tstamp
            row = pd.DataFrame(row)
            data = pd.concat([data, row])

    cal_dicts, results_dicts, object_dicts, plot_dicts = calibrate_partition(
        data,
        cal_dicts,
        results_dicts,
        object_dicts,
        inplots_dict,
        args.channel,
        chmap,
        kwarg_dict,
        gen_plots=bool(args.plot_file),
        debug_mode=args.debug,
    )

    if args.plot_file:
        save_dict_to_files(args.plot_file, plot_dicts)

    save_dict_to_files(
        args.hit_pars,
        {
            tstamp: {"pars": cal_dicts[tstamp], "results": results_dicts[tstamp]}
            for tstamp in cal_dicts
        },
    )
    save_dict_to_files(args.fit_results, object_dicts)
