from __future__ import annotations

import argparse
import warnings

import numpy as np
import pandas as pd
from dbetto import Props, TextDB
from legenddataflowscripts.utils import build_log, get_pulser_mask
from pygama.pargen.AoE_cal import *  # noqa: F403
from pygama.pargen.lq_cal import *  # noqa: F403
from pygama.pargen.utils import load_data

from legenddataflow.scripts.par.geds.hit.lq import run_lq_calibration

warnings.filterwarnings(action="ignore", category=RuntimeWarning)


def par_geds_pht_lq() -> None:
    argparser = argparse.ArgumentParser()
    argparser.add_argument(
        "--input-files", help="files", type=str, nargs="*", required=True
    )
    argparser.add_argument(
        "--pulser-files", help="pulser_file", type=str, nargs="*", required=False
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
    argparser.add_argument("--metadata", help="metadata path", type=str, required=True)
    argparser.add_argument("--log", help="log_file", type=str)

    argparser.add_argument("--timestamp", help="Datatype", type=str, required=True)
    argparser.add_argument("--datatype", help="Datatype", type=str, required=True)
    argparser.add_argument("--channel", help="Channel", type=str, required=True)
    argparser.add_argument("--table-name", help="table name", type=str, required=True)

    argparser.add_argument(
        "--plot-file", help="plot_file", type=str, nargs="*", required=False
    )
    argparser.add_argument("--hit-pars", help="hit_pars", nargs="*", type=str)
    argparser.add_argument("--lq-results", help="lq_results", nargs="*", type=str)

    argparser.add_argument("-d", "--debug", help="debug_mode", action="store_true")
    args = argparser.parse_args()

    configs = TextDB(args.configs, lazy=True).on(args.timestamp, system=args.datatype)
    config_dict = configs["snakemake_rules"]["pars_pht_lqcal"]

    log = build_log(config_dict, args.log)

    channel_dict = config_dict["inputs"]["lqcal_config"][args.channel]
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
            args.table_name,
            cal_dicts,
            params=params,
            threshold=kwarg_dict.pop("threshold"),
            return_selection_mask=True,
        )
        msg = f"Loaded {len(data)} events from {len(final_dict)} files"
        log.info(msg)

        mask = get_pulser_mask(pulser_file=args.pulser_files)
        if "pulser_multiplicity_threshold" in kwarg_dict:
            kwarg_dict.pop("pulser_multiplicity_threshold")

        data["is_pulser"] = mask[threshold_mask]
        msg = f"{len(data.query('~is_pulser'))}  non pulser events"
        log.info(msg)

        for tstamp in cal_dicts:
            if tstamp not in np.unique(data["run_timestamp"]):
                row = {
                    key: [False] if data.dtypes[key] == "bool" else [np.nan]
                    for key in data
                }
                row["run_timestamp"] = tstamp
                row = pd.DataFrame(row)
                data = pd.concat([data, row])

        cal_dicts, results_dicts, object_dicts, plot_dicts = run_lq_calibration(
            data,
            cal_dicts,
            results_dicts,
            object_dicts,
            inplots_dict,
            kwarg_dict,
            # gen_plots=bool(args.plot_file),
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
    save_dict_to_files(args.lq_results, object_dicts)
