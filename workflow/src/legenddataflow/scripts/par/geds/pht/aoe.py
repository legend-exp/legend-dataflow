from __future__ import annotations

import argparse
import copy
import pickle as pkl
import re
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from dbetto import TextDB
from dbetto.catalog import Props
from legendmeta import LegendMetadata
from pygama.pargen.AoE_cal import *  # noqa: F403
from pygama.pargen.AoE_cal import CalAoE, Pol1, SigmaFit, aoe_peak
from pygama.pargen.utils import load_data

from .....FileKey import ChannelProcKey, ProcessingFileKey
from .....log import build_log
from ....pulser_removal import get_pulser_mask

warnings.filterwarnings(action="ignore", category=RuntimeWarning)


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
            "low_side_sfs_by_run": aoe_class.low_side_sfs_by_run[tstamp].to_dict(
                "index"
            ),
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


def run_aoe_calibration(
    data,
    cal_dicts,
    results_dicts,
    object_dicts,
    plot_dicts,
    timestamp,
    configs,
    channel,
    datatype,
    debug_mode=False,
    # gen_plots=True,
):
    configs = LegendMetadata(path=configs)
    channel_dict = configs.on(timestamp, system=datatype)["snakemake_rules"][
        "pars_pht_aoecal"
    ]["inputs"]["par_pht_aoecal_config"][channel]

    kwarg_dict = Props.read_from(channel_dict)

    if kwarg_dict.pop("run_aoe") is True:
        kwarg_dict.pop("pulser_multiplicity_threshold")
        kwarg_dict.pop("threshold")
        if "plot_options" in kwarg_dict:
            for field, item in kwarg_dict["plot_options"].items():
                kwarg_dict["plot_options"][field]["function"] = eval(item["function"])

        pdf = eval(kwarg_dict.pop("pdf")) if "pdf" in kwarg_dict else aoe_peak

        mean_func = (
            eval(kwarg_dict.pop("mean_func")) if "mean_func" in kwarg_dict else Pol1
        )

        sigma_func = (
            eval(kwarg_dict.pop("sigma_func"))
            if "sigma_func" in kwarg_dict
            else SigmaFit
        )

        if "dt_cut" in kwarg_dict and kwarg_dict["dt_cut"] is not None:
            cut_dict = kwarg_dict["dt_cut"]["cut"]
            for tstamp in cal_dicts:
                cal_dicts[tstamp].update(cut_dict)

            exp = cut_dict[next(iter(cut_dict))]["expression"]
            for key in cut_dict[next(iter(cut_dict))]["parameters"]:
                exp = re.sub(f"(?<![a-zA-Z0-9]){key}(?![a-zA-Z0-9])", f"@{key}", exp)
            data[next(iter(cut_dict))] = data.eval(
                exp, local_dict=cut_dict[next(iter(cut_dict))]["parameters"]
            )

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
                    results_dicts[next(iter(results_dicts))]["ecal"][
                        kwarg_dict["cal_energy_param"]
                    ]["eres_linear"]
                )

                def eres_func(x):
                    return eval(eres["expression"], dict(x=x, **eres["parameters"]))

            except KeyError:

                def eres_func(x):
                    return x * np.nan

        data["AoE_Uncorr"] = (
            data[kwarg_dict["current_param"]] / data[kwarg_dict["energy_param"]]
        )
        aoe = CalAoE(
            cal_dicts=cal_dicts,
            cal_energy_param=kwarg_dict["cal_energy_param"],
            eres_func=eres_func,
            pdf=pdf,
            mean_func=mean_func,
            sigma_func=sigma_func,
            selection_string=f"{kwarg_dict.pop('final_cut_field')}&(~is_pulser)",
            dt_corr=kwarg_dict.get("dt_corr", False),
            dep_correct=kwarg_dict.get("dep_correct", False),
            dt_cut=kwarg_dict.get("dt_cut", None),
            dt_param=kwarg_dict.get("dt_param", 3),
            high_cut_val=kwarg_dict.get("high_cut_val", 3),
            compt_bands_width=kwarg_dict.get("debug_mode", 20),
            debug_mode=debug_mode | kwarg_dict.get("debug_mode", False),
        )
        aoe.update_cal_dicts(
            {
                "AoE_Uncorr": {
                    "expression": f"{kwarg_dict['current_param']}/{kwarg_dict['energy_param']}",
                    "parameters": {},
                }
            }
        )
        aoe.calibrate(data, "AoE_Uncorr")

        out_dict = get_results_dict(aoe)
        aoe_plot_dict = fill_plot_dict(aoe, data, kwarg_dict.get("plot_options", None))

        aoe.pdf = aoe.pdf.name
        # need to change eres func as can't pickle lambdas
        try:
            aoe.eres_func = results_dicts[next(iter(results_dicts))]["partition_ecal"][
                kwarg_dict["cal_energy_param"]
            ]["eres_linear"]
        except KeyError:
            aoe.eres_func = {}
    else:
        out_dict = {tstamp: None for tstamp in cal_dicts}
        aoe_plot_dict = {}
        aoe = None

    out_result_dicts = {}
    for tstamp, result_dict in results_dicts.items():
        out_result_dicts[tstamp] = dict(**result_dict, aoe=out_dict[tstamp])

    out_object_dicts = {}
    for tstamp, object_dict in object_dicts.items():
        out_object_dicts[tstamp] = dict(**object_dict, aoe=aoe)

    common_dict = (
        aoe_plot_dict.pop("common") if "common" in list(aoe_plot_dict) else None
    )
    out_plot_dicts = {}
    for tstamp, plot_dict in plot_dicts.items():
        if "common" in list(plot_dict) and common_dict is not None:
            plot_dict["common"].update(common_dict)
        elif common_dict is not None:
            plot_dict["common"] = common_dict
        plot_dict.update({"aoe": aoe_plot_dict})
        out_plot_dicts[tstamp] = plot_dict

    return cal_dicts, out_result_dicts, out_object_dicts, out_plot_dicts


def par_geds_pht_aoe() -> None:
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

    argparser.add_argument("--configs", help="configs", type=str, required=True)
    argparser.add_argument("--metadata", help="metadata", type=str)
    argparser.add_argument("--log", help="log_file", type=str)

    argparser.add_argument("--timestamp", help="Datatype", type=str, required=True)
    argparser.add_argument("--datatype", help="Datatype", type=str, required=True)
    argparser.add_argument("--channel", help="Channel", type=str, required=True)
    argparser.add_argument("--table-name", help="table name", type=str, required=True)

    argparser.add_argument(
        "--plot-file", help="plot_file", type=str, nargs="*", required=False
    )
    argparser.add_argument("--hit-pars", help="hit_pars", nargs="*", type=str)
    argparser.add_argument("--aoe-results", help="aoe_results", nargs="*", type=str)

    argparser.add_argument("-d", "--debug", help="debug_mode", action="store_true")
    args = argparser.parse_args()

    configs = TextDB(args.configs, lazy=True).on(args.timestamp, system=args.datatype)
    config_dict = configs["snakemake_rules"]["pars_pht_aoecal"]

    build_log(config_dict, args.log)

    channel_dict = config_dict["inputs"]["par_pht_aoecal_config"][args.channel]
    kwarg_dict = Props.read_from(channel_dict)

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
    if isinstance(args.input_files, list):
        files = []
        for file in args.input_files:
            with Path(file).open() as f:
                files += f.read().splitlines()
    else:
        with Path(args.input_files).open() as f:
            files = f.read().splitlines()

    files = sorted(
        np.unique(files)
    )  # need this as sometimes files get double counted as it somehow puts in the p%-* filelist and individual runs also

    final_dict = {}
    all_file = run_splitter(sorted(files))
    for filelist in all_file:
        fk = ProcessingFileKey.get_filekey_from_pattern(Path(sorted(filelist)[0]).name)
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

        # load data in
        data, threshold_mask = load_data(
            final_dict,
            args.table_name,
            cal_dict,
            params=params,
            threshold=kwarg_dict.pop("threshold"),
            return_selection_mask=True,
        )

        mask = get_pulser_mask(pulser_file=args.pulser_files)
        if "pulser_multiplicity_threshold" in kwarg_dict:
            kwarg_dict.pop("pulser_multiplicity_threshold")

        data["is_pulser"] = mask[threshold_mask]

        for tstamp in cal_dict:
            if tstamp not in np.unique(data["run_timestamp"]):
                row = {
                    key: [False] if data.dtypes[key] == "bool" else [np.nan]
                    for key in data
                }
                row["run_timestamp"] = tstamp
                row = pd.DataFrame(row)
                data = pd.concat([data, row])

        cal_dicts, results_dicts, object_dicts, plot_dicts = run_aoe_calibration(
            data,
            cal_dict,
            results_dicts,
            object_dict,
            inplots_dict,
            timestamp,
            args.configs,
            args.channel,
            args.datatype,
            debug_mode=args.debug,
            # gen_plots=bool(args.plot_file),
        )

        if args.plot_file:
            for plot_file in args.plot_file:
                Path(plot_file).parent.mkdir(parents=True, exist_ok=True)
                with Path(plot_file).open("wb") as w:
                    pkl.dump(plot_dicts[fk.timestamp], w, protocol=pkl.HIGHEST_PROTOCOL)

        for out in sorted(args.hit_pars):
            fk = ChannelProcKey.get_filekey_from_pattern(Path(out).name)
            final_hit_dict = {
                "pars": cal_dicts[fk.timestamp],
                "results": results_dicts[fk.timestamp],
            }
            Path(out).parent.mkdir(parents=True, exist_ok=True)
            Props.write_to(out, final_hit_dict)

        for out in args.aoe_results:
            fk = ChannelProcKey.get_filekey_from_pattern(Path(out).name)
            Path(out).parent.mkdir(parents=True, exist_ok=True)
            with Path(out).open("wb") as w:
                pkl.dump(object_dicts[fk.timestamp], w, protocol=pkl.HIGHEST_PROTOCOL)
