from __future__ import annotations

import argparse
import warnings

import numpy as np
import pandas as pd
from dbetto import Props, TextDB
from legenddataflow.scripts.par.geds.hit.aoe import run_aoe_calibration
from pygama.pargen.AoE_cal import *  # noqa: F403
from pygama.pargen.utils import load_data

from .....log import build_log
from ....pulser_removal import get_pulser_mask
from .util import get_run_dict, save_dict_to_files, split_files_by_run

warnings.filterwarnings(action="ignore", category=RuntimeWarning)


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

    log = build_log(config_dict, args.log)

    channel_dict = config_dict["inputs"]["par_pht_aoecal_config"][args.channel]
    kwarg_dict = Props.read_from(channel_dict)

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

    # run aoe cal
    if kwarg_dict.get("run_aoe", True) is True:
        params = [
            kwarg_dict["final_cut_field"],
            kwarg_dict["current_param"],
            "tp_0_est",
            "tp_99",
            kwarg_dict["energy_param"],
            kwarg_dict["cal_energy_param"],
            "timestamp",
        ]
        params.append(kwarg_dict.get("dt_param", "dt_eff"))

        if "dt_cut" in kwarg_dict and kwarg_dict["dt_cut"] is not None:
            params.append(kwarg_dict["dt_cut"]["out_param"])

        # load data in
        data, threshold_mask = load_data(
            final_dict,
            args.table_name,
            cal_dicts,
            params=params,
            threshold=kwarg_dict.get("threshold", 900),
            return_selection_mask=True,
        )

        mask = get_pulser_mask(pulser_file=args.pulser_files)
        data["is_pulser"] = mask[threshold_mask]

        for tstamp in cal_dicts:
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
            cal_dicts,
            results_dicts,
            object_dicts,
            configs=channel_dict,
            log=log,
            debug_mode=args.debug,
            # gen_plots=bool(args.plot_file),
        )
    else:
        plot_dicts = inplots_dict

    save_dict_to_files(
        args.hit_pars,
        {
            tstamp: {"pars": cal_dicts[tstamp], "results": results_dicts[tstamp]}
            for tstamp in cal_dicts
        },
    )
    save_dict_to_files(args.aoe_results, object_dicts)
    if args.plot_file:
        save_dict_to_files(args.plot_file, plot_dicts)
