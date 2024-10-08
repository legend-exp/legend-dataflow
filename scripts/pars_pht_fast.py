from __future__ import annotations

import argparse
import json
import logging
import os
import pathlib
import pickle as pkl
import warnings

import numpy as np
import pandas as pd
from legendmeta import LegendMetadata
from legendmeta.catalog import Props
from pars_pht_aoecal import run_aoe_calibration
from pars_pht_lqcal import run_lq_calibration
from pars_pht_partcal import calibrate_partition
from pygama.pargen.data_cleaning import get_tcm_pulser_ids
from pygama.pargen.utils import load_data
from util.FileKey import ChannelProcKey, ProcessingFileKey

log = logging.getLogger(__name__)
warnings.filterwarnings(action="ignore", category=RuntimeWarning)
warnings.filterwarnings(action="ignore", category=np.RankWarning)


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

    cal_dict = {}
    results_dicts = {}
    for ecal in args.ecal_file:
        cal = Props.read_from(ecal)

        fk = ChannelProcKey.get_filekey_from_pattern(os.path.basename(ecal))
        cal_dict[fk.timestamp] = cal["pars"]
        results_dicts[fk.timestamp] = cal["results"]

    object_dict = {}
    for ecal in args.eres_file:
        with open(ecal, "rb") as o:
            cal = pkl.load(o)
        fk = ChannelProcKey.get_filekey_from_pattern(os.path.basename(ecal))
        object_dict[fk.timestamp] = cal

    inplots_dict = {}
    if args.inplots:
        for ecal in args.inplots:
            with open(ecal, "rb") as o:
                cal = pkl.load(o)
            fk = ChannelProcKey.get_filekey_from_pattern(os.path.basename(ecal))
            inplots_dict[fk.timestamp] = cal

    # sort files in dictionary where keys are first timestamp from run
    files = []
    for file in args.input_files:
        with open(file) as f:
            files += f.read().splitlines()

    files = sorted(
        np.unique(files)
    )  # need this as sometimes files get double counted as it somehow puts in the p%-* filelist and individual runs also

    final_dict = {}
    all_file = run_splitter(sorted(files))
    for filelist in all_file:
        fk = ProcessingFileKey.get_filekey_from_pattern(os.path.basename(sorted(filelist)[0]))
        timestamp = fk.timestamp
        final_dict[timestamp] = sorted(filelist)

    configs = LegendMetadata(path=args.configs)
    channel_dict = configs.on(timestamp, system=args.datatype)["snakemake_rules"]

    kwarg_dict = Props.read_from(
        channel_dict["pars_pht_partcal"]["inputs"]["pars_pht_partcal_config"][args.channel]
    )
    aoe_kwarg_dict = Props.read_from(
        channel_dict["pars_pht_aoecal"]["inputs"]["par_pht_aoecal_config"][args.channel]
    )
    lq_kwarg_dict = Props.read_from(
        channel_dict["pars_pht_lqcal"]["inputs"]["lqcal_config"][args.channel]
    )

    params = [
        kwarg_dict["final_cut_field"],
        "timestamp",
    ]
    params += kwarg_dict["energy_params"]

    if aoe_kwarg_dict["run_aoe"] is True:
        aoe_params = [
            aoe_kwarg_dict["final_cut_field"],
            aoe_kwarg_dict["current_param"],
            "tp_0_est",
            "tp_99",
            aoe_kwarg_dict["energy_param"],
            aoe_kwarg_dict["cal_energy_param"],
            "timestamp",
        ]
        if "dt_param" in aoe_kwarg_dict:
            aoe_params.append(aoe_kwarg_dict["dt_param"])
        else:
            aoe_params.append("dt_eff")

        params += aoe_params

    if lq_kwarg_dict["run_lq"] is True:
        params += [
            "lq80",
            "dt_eff",
            lq_kwarg_dict["energy_param"],
            lq_kwarg_dict["cal_energy_param"],
            lq_kwarg_dict["cut_field"],
        ]
    params = list(np.unique(params))

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

    cal_dicts, results_dicts, object_dicts, plot_dicts = calibrate_partition(
        data,
        cal_dict,
        results_dicts,
        object_dict,
        inplots_dict,
        args.timestamp,
        args.metadata,
        args.configs,
        args.channel,
        args.datatype,
        gen_plots=bool(args.plot_file),
    )

    cal_dicts, results_dicts, object_dicts, plot_dicts = run_aoe_calibration(
        data,
        cal_dicts,
        results_dicts,
        object_dicts,
        plot_dicts,
        args.timestamp,
        args.configs,
        args.channel,
        args.datatype,
        # gen_plots=bool(args.plot_file),
    )

    cal_dicts, results_dicts, object_dicts, plot_dicts = run_lq_calibration(
        data,
        cal_dicts,
        results_dicts,
        object_dicts,
        plot_dicts,
        args.timestamp,
        args.configs,
        args.channel,
        args.datatype,
        # gen_plots=bool(args.plot_file),
    )

    if args.plot_file:
        for plot_file in args.plot_file:
            pathlib.Path(os.path.dirname(plot_file)).mkdir(parents=True, exist_ok=True)
            with open(plot_file, "wb") as w:
                pkl.dump(plot_dicts[fk.timestamp], w, protocol=pkl.HIGHEST_PROTOCOL)

    for out in sorted(args.hit_pars):
        fk = ChannelProcKey.get_filekey_from_pattern(os.path.basename(out))
        final_hit_dict = {
            "pars": {"operations": cal_dict[fk.timestamp]},
            "results": results_dicts[fk.timestamp],
        }
        pathlib.Path(os.path.dirname(out)).mkdir(parents=True, exist_ok=True)
        with open(out, "w") as w:
            json.dump(final_hit_dict, w, indent=4)

    for out in args.fit_results:
        fk = ChannelProcKey.get_filekey_from_pattern(os.path.basename(out))
        pathlib.Path(os.path.dirname(out)).mkdir(parents=True, exist_ok=True)
        with open(out, "wb") as w:
            pkl.dump(object_dicts[fk.timestamp], w, protocol=pkl.HIGHEST_PROTOCOL)
