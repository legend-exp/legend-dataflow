from __future__ import annotations

import argparse
import copy
import json
import pickle as pkl
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from dbetto import TextDB
from dbetto.catalog import Props
from legendmeta import LegendMetadata
from pygama.math.distributions import gaussian
from pygama.pargen.AoE_cal import *  # noqa: F403
from pygama.pargen.data_cleaning import get_tcm_pulser_ids
from pygama.pargen.lq_cal import *  # noqa: F403
from pygama.pargen.lq_cal import LQCal
from pygama.pargen.utils import load_data

from ..FileKey import ChannelProcKey, ProcessingFileKey
from ..log import build_log

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


def get_results_dict(lq_class):
    return {
        "cal_energy_param": lq_class.cal_energy_param,
        "DEP_means": lq_class.timecorr_df.to_dict("index"),
        "rt_correction": lq_class.dt_fit_pars,
        "cut_fit_pars": lq_class.cut_fit_pars.to_dict(),
        "cut_value": lq_class.cut_val,
        "sfs": lq_class.low_side_sf.to_dict("index"),
    }


def fill_plot_dict(lq_class, data, plot_options, plot_dict=None):
    if plot_dict is None:
        plot_dict = {}
    for key, item in plot_options.items():
        if item["options"] is not None:
            plot_dict[key] = item["function"](lq_class, data, **item["options"])
        else:
            plot_dict[key] = item["function"](lq_class, data)

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
        The drift-time parameter of choice
    eres_func: callable
        The energy resolution functions
    cdf: callable
        The CDF used for the binned fitting of LQ distributions
    cut_field: string
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
    lq: cal_lq class
        The cal_lq object used for the LQ calibration
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


def run_lq_calibration(
    data,
    cal_dicts,
    results_dicts,
    object_dicts,
    plot_dicts,
    timestamp,
    configs,
    channel,
    datatype,
    # gen_plots=True,
):
    configs = LegendMetadata(path=configs)
    channel_dict = configs.on(timestamp, system=datatype)["snakemake_rules"]["pars_pht_lqcal"][
        "inputs"
    ]["lqcal_config"][channel]

    kwarg_dict = Props.read_from(channel_dict)

    if kwarg_dict.pop("run_lq") is True:

        if "plot_options" in kwarg_dict:
            for field, item in kwarg_dict["plot_options"].items():
                kwarg_dict["plot_options"][field]["function"] = eval(item["function"])

        kwarg_dict.pop("pulser_multiplicity_threshold")
        kwarg_dict.pop("threshold")

        cdf = eval(kwarg_dict.pop("cdf")) if "cdf" in kwarg_dict else gaussian

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

        cal_dicts, out_dict, lq_plot_dict, lq_obj = lq_calibration(
            data,
            selection_string=f"{kwarg_dict.pop('cut_field')}&(~is_pulser)",
            cal_dicts=cal_dicts,
            eres_func=eres_func,
            cdf=cdf,
            **kwarg_dict,
        )
        # need to change eres func as can't pickle lambdas
        try:
            lq_obj.eres_func = results_dicts[next(iter(results_dicts))]["partition_ecal"][
                kwarg_dict["cal_energy_param"]
            ]["eres_linear"]
        except KeyError:
            lq_obj.eres_func = {}
    else:
        out_dict = {tstamp: None for tstamp in cal_dicts}
        lq_plot_dict = {}
        lq_obj = None

    out_result_dicts = {}
    for tstamp, result_dict in results_dicts.items():
        out_result_dicts[tstamp] = dict(**result_dict, lq=out_dict)

    out_object_dicts = {}
    for tstamp, object_dict in object_dicts.items():
        out_object_dicts[tstamp] = dict(**object_dict, lq=lq_obj)

    common_dict = lq_plot_dict.pop("common") if "common" in list(lq_plot_dict) else None
    out_plot_dicts = {}
    for tstamp, plot_dict in plot_dicts.items():
        if "common" in list(plot_dict) and common_dict is not None:
            plot_dict["common"].update(lq_plot_dict["common"])
        elif common_dict is not None:
            plot_dict["common"] = common_dict
        plot_dict.update({"lq": lq_plot_dict})
        out_plot_dicts[tstamp] = plot_dict

    return cal_dicts, out_result_dicts, out_object_dicts, out_plot_dicts


if __name__ == "__main__":
    argparser = argparse.ArgumentParser()
    argparser.add_argument("--input_files", help="files", type=str, nargs="*", required=True)
    argparser.add_argument(
        "--pulser_files", help="pulser_file", type=str, nargs="*", required=False
    )
    argparser.add_argument(
        "--tcm_filelist", help="tcm_filelist", type=str, nargs="*", required=False
    )
    argparser.add_argument("--ecal_file", help="ecal_file", type=str, nargs="*", required=True)
    argparser.add_argument("--eres_file", help="eres_file", type=str, nargs="*", required=True)
    argparser.add_argument("--inplots", help="eres_file", type=str, nargs="*", required=True)

    argparser.add_argument("--configs", help="configs", type=str, required=True)
    argparser.add_argument("--metadata", help="metadata path", type=str, required=True)
    argparser.add_argument("--log", help="log_file", type=str)

    argparser.add_argument("--timestamp", help="Datatype", type=str, required=True)
    argparser.add_argument("--datatype", help="Datatype", type=str, required=True)
    argparser.add_argument("--channel", help="Channel", type=str, required=True)

    argparser.add_argument("--plot_file", help="plot_file", type=str, nargs="*", required=False)
    argparser.add_argument("--hit_pars", help="hit_pars", nargs="*", type=str)
    argparser.add_argument("--lq_results", help="lq_results", nargs="*", type=str)

    argparser.add_argument("-d", "--debug", help="debug_mode", action="store_true")
    args = argparser.parse_args()

    configs = TextDB(args.configs, lazy=True).on(args.timestamp, system=args.datatype)
    config_dict = configs["snakemake_rules"]["pars_pht_lqcal"]

    log = build_log(config_dict, args.log)

    meta = LegendMetadata(path=args.metadata)
    chmap = meta.channelmap(args.timestamp, system=args.datatype)
    channel = f"ch{chmap[args.channel].daq.rawid:07}"

    channel_dict = config_dict["inputs"]["lqcal_config"][args.channel]
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

    # run lq cal
    if kwarg_dict.pop("run_lq") is True:
        params = [
            "lq80",
            "dt_eff",
            kwarg_dict["energy_param"],
            kwarg_dict["cal_energy_param"],
            kwarg_dict["cut_field"],
        ]

        # load data in
        data, threshold_mask = load_data(
            final_dict,
            f"{channel}/dsp",
            cal_dict,
            params=params,
            threshold=kwarg_dict.pop("threshold"),
            return_selection_mask=True,
        )

        if args.pulser_files:
            mask = np.array([], dtype=bool)
            for file in args.pulser_files:
                with Path(file).open() as f:
                    pulser_dict = json.load(f)
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

        cal_dicts, results_dicts, object_dicts, plot_dicts = run_lq_calibration(
            data,
            cal_dict,
            results_dicts,
            object_dict,
            inplots_dict,
            timestamp,
            args.configs,
            args.channel,
            args.datatype,
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
            "pars": {"operations": cal_dict[fk.timestamp]},
            "results": results_dicts[fk.timestamp],
        }
        Path(out).parent.mkdir(parents=True, exist_ok=True)
        with Path(out).open("w") as w:
            json.dump(final_hit_dict, w, indent=4)

    for out in args.lq_results:
        fk = ChannelProcKey.get_filekey_from_pattern(Path(out).name)
        Path(out).parent.mkdir(parents=True, exist_ok=True)
        with Path(out).open("wb") as w:
            pkl.dump(object_dict[fk.timestamp], w, protocol=pkl.HIGHEST_PROTOCOL)
