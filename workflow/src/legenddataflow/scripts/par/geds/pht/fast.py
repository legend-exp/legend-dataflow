from __future__ import annotations

import argparse
import json
import pickle as pkl
import time
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from dbetto import TextDB
from dbetto.catalog import Props
from legenddataflowscripts.utils import build_log, get_pulser_mask
from legendmeta import LegendMetadata
from pygama.pargen.utils import load_data

from legenddataflow.methods import ChannelProcKey, ProcessingFileKey, run_splitter

from .aoe import run_aoe_calibration
from .ecal_part import calibrate_partition
from .lq import run_lq_calibration

warnings.filterwarnings(action="ignore", category=RuntimeWarning)
warnings.filterwarnings(action="ignore", category=np.exceptions.RankWarning)


def par_geds_pht_fast() -> None:
    argparser = argparse.ArgumentParser()
    argparser.add_argument(
        "--input-files", help="files", type=str, nargs="*", required=True
    )
    argparser.add_argument(
        "--pulser-files", help="pulser_file", nargs="*", type=str, required=False
    )
    argparser.add_argument(
        "--tcm-filelist", help="tcm_filelist", type=str, nargs="*", required=False
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
    config_dict = configs["snakemake_rules"]

    log = build_log(config_dict["pars_pht_partcal"], args.log)

    chmap = LegendMetadata(args.metadata).channelmap(
        args.timestamp, system=args.datatype
    )

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

    kwarg_dict = Props.read_from(
        config_dict["pars_pht_partcal"]["inputs"]["pars_pht_partcal_config"][
            args.channel
        ]
    )
    aoe_kwarg_dict = Props.read_from(
        config_dict["pars_pht_aoecal"]["inputs"]["par_pht_aoecal_config"][args.channel]
    )
    lq_kwarg_dict = Props.read_from(
        config_dict["pars_pht_lqcal"]["inputs"]["lqcal_config"][args.channel]
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
        args.table_name,
        cal_dict,
        params=params,
        threshold=kwarg_dict["threshold"],
        return_selection_mask=True,
        cal_energy_param=kwarg_dict["energy_params"][0],
    )
    msg = f"Loaded {len(data)} events"
    log.info(msg)

    mask = get_pulser_mask(pulser_file=args.pulser_files)
    if "pulser_multiplicity_threshold" in kwarg_dict:
        kwarg_dict.pop("pulser_multiplicity_threshold")

    data["is_pulser"] = mask[threshold_mask]

    msg = f"{len(data.query('~is_pulser'))} non pulser events"
    log.info(msg)

    for tstamp in cal_dict:
        if tstamp not in np.unique(data["run_timestamp"]):
            row = {
                key: [False] if data.dtypes[key] == "bool" else [np.nan] for key in data
            }
            row["run_timestamp"] = tstamp
            row = pd.DataFrame(row)
            data = pd.concat([data, row])

    configs = TextDB(path=args.configs, lazy=True).on(
        args.timestamp, system=args.datatype
    )["snakemake_rules"]

    start = time.time()

    cal_dicts, results_dicts, object_dicts, plot_dicts = calibrate_partition(
        data,
        cal_dict,
        results_dicts,
        object_dict,
        inplots_dict,
        args.timestamp,
        chmap,
        config=configs["pars_pht_aoecal"]["inputs"]["par_pht_partcal_config"][
            args.channel
        ],
        gen_plots=bool(args.plot_file),
        debug_mode=args.debug,
    )
    start2 = time.time()
    msg = f"Partition calibration took {start2 - start:.2f} seconds"
    log.info(msg)

    cal_dicts, results_dicts, object_dicts, plot_dicts = run_aoe_calibration(
        data,
        cal_dicts,
        results_dicts,
        object_dicts,
        plot_dicts,
        config=configs["pars_pht_aoecal"]["inputs"]["par_pht_aoecal_config"][
            args.channel
        ],
        debug_mode=args.debug,
        # gen_plots=bool(args.plot_file),
    )

    start3 = time.time()
    msg = f"A/E calibration took {start3 - start2:.2f} seconds"
    log.info(msg)

    cal_dicts, results_dicts, object_dicts, plot_dicts = run_lq_calibration(
        data,
        cal_dicts,
        results_dicts,
        object_dicts,
        plot_dicts,
        config=configs["pars_pht_aoecal"]["inputs"]["par_pht_lqcal_config"][
            args.channel
        ],
        debug_mode=args.debug,
        # gen_plots=bool(args.plot_file),
    )

    msg = f"LQ calibration took {time.time() - start3:.2f} seconds"
    log.info(msg)
    msg = f"Total calibration took {time.time() - start:.2f} seconds"
    log.info(msg)

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

    for out in args.fit_results:
        fk = ChannelProcKey.get_filekey_from_pattern(Path(out).name)
        Path(out).parent.mkdir(parents=True, exist_ok=True)
        with Path(out).open("wb") as w:
            pkl.dump(object_dicts[fk.timestamp], w, protocol=pkl.HIGHEST_PROTOCOL)
