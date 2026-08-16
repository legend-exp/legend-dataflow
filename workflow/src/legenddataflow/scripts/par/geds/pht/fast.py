"""Console script ``par-geds-pht-fast``: run the partition-level calibrations
(energy, A/E, LQ) for a HPGe channel in a single pass."""

from __future__ import annotations

import argparse
import pickle as pkl
import time
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from dbetto.catalog import Props
from legenddataflowscripts.par.geds.hit.aoe import run_aoe_calibration
from legenddataflowscripts.par.geds.hit.lq import run_lq_calibration
from legenddataflowscripts.utils import (
    build_log,
    check_input_files,
    check_pulser_mask,
    get_channel_config,
    get_pulser_mask,
    get_rule_config,
    prepare_output_paths,
    require_config_keys,
)
from legendmeta import LegendMetadata
from pygama.pargen.utils import load_data

from legenddataflow.methods import ChannelProcKey, ProcessingFileKey, run_splitter

from .ecal_part import calibrate_partition
from .util import (
    save_dict_to_files,
)

warnings.filterwarnings(action="ignore", category=RuntimeWarning)
try:
    warnings.filterwarnings(action="ignore", category=np.exceptions.RankWarning)
except AttributeError:  # np < 2
    warnings.filterwarnings(action="ignore", category=np.RankWarning)


def par_geds_pht_fast() -> None:
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

    partcal_dict = get_rule_config(
        args.configs, "pars_pht_partcal", args.timestamp, args.datatype
    )
    aoecal_dict = get_rule_config(
        args.configs, "pars_pht_aoecal", args.timestamp, args.datatype
    )
    lqcal_dict = get_rule_config(
        args.configs, "pars_pht_lqcal", args.timestamp, args.datatype
    )

    log = build_log(partcal_dict, args.log)

    check_input_files(args.pulser_files, "--pulser-files")
    prepare_output_paths(
        *(args.hit_pars or []), *(args.fit_results or []), *(args.plot_file or [])
    )

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
    check_input_files(files, "--input-files")

    final_dict = {}
    all_file = run_splitter(sorted(files))
    for filelist in all_file:
        fk = ProcessingFileKey.get_filekey_from_pattern(Path(sorted(filelist)[0]).name)
        timestamp = fk.timestamp
        final_dict[timestamp] = sorted(filelist)

    partcal_config = get_channel_config(
        partcal_dict["inputs"]["pars_pht_partcal_config"],
        args.channel,
        name="pars_pht_partcal_config",
    )
    aoecal_config = get_channel_config(
        aoecal_dict["inputs"]["par_pht_aoecal_config"],
        args.channel,
        name="par_pht_aoecal_config",
    )
    lqcal_config = get_channel_config(
        lqcal_dict["inputs"]["lqcal_config"],
        args.channel,
        name="lqcal_config",
    )
    kwarg_dict = Props.read_from(partcal_config)
    aoe_kwarg_dict = Props.read_from(aoecal_config)
    lq_kwarg_dict = Props.read_from(lqcal_config)

    require_config_keys(
        kwarg_dict,
        ["final_cut_field", "energy_params", "threshold"],
        f"channel {args.channel} partcal config ({partcal_config})",
    )
    require_config_keys(
        aoe_kwarg_dict,
        ["run_aoe"],
        f"channel {args.channel} aoecal config ({aoecal_config})",
    )
    require_config_keys(
        lq_kwarg_dict,
        ["run_lq"],
        f"channel {args.channel} lqcal config ({lqcal_config})",
    )

    params = [
        kwarg_dict["final_cut_field"],
        "timestamp",
    ]
    params += kwarg_dict["energy_params"]

    if aoe_kwarg_dict["run_aoe"] is True:
        require_config_keys(
            aoe_kwarg_dict,
            ["final_cut_field", "cal_energy_param", "cut_field", "params"],
            f"channel {args.channel} aoecal config ({aoecal_config})",
        )
        aoe_params = [
            aoe_kwarg_dict["final_cut_field"],
            aoe_kwarg_dict["cal_energy_param"],
            aoe_kwarg_dict["cut_field"],
            "tp_0_est",
            "tp_99",
            "timestamp",
        ]
        for name, param_config in aoe_kwarg_dict["params"].items():
            require_config_keys(
                param_config,
                ["current_param", "energy_param"],
                f"channel {args.channel} aoecal config params entry '{name}' ({aoecal_config})",
            )
            aoe_params.append(param_config["current_param"])
            aoe_params.append(param_config["energy_param"])
            aoe_params.append(param_config.get("dt_param", "dt_eff"))
            if "dt_cut" in param_config and param_config["dt_cut"] is not None:
                require_config_keys(
                    param_config["dt_cut"],
                    ["out_param"],
                    f"channel {args.channel} aoecal config params entry '{name}' dt_cut ({aoecal_config})",
                )
                aoe_params.append(param_config["dt_cut"]["out_param"])

        params += aoe_params

    if lq_kwarg_dict["run_lq"] is True:
        require_config_keys(
            lq_kwarg_dict,
            ["cal_energy_param", "cut_field", "params"],
            f"channel {args.channel} lqcal config ({lqcal_config})",
        )
        params += [
            lq_kwarg_dict["cal_energy_param"],
            lq_kwarg_dict["cut_field"],
        ]
        for name, param_config in lq_kwarg_dict["params"].items():
            require_config_keys(
                param_config,
                ["lq_param", "energy_param"],
                f"channel {args.channel} lqcal config params entry '{name}' ({lqcal_config})",
            )
            params.append(param_config["lq_param"])
            params.append(param_config["energy_param"])
            params.append(param_config.get("dt_param", "dt_eff"))
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

    check_pulser_mask(mask, threshold_mask, args.channel)
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

    start = time.time()

    cal_dicts, results_dicts, object_dicts, plot_dicts = calibrate_partition(
        data,
        cal_dict,
        results_dicts,
        object_dict,
        inplots_dict,
        args.channel,
        chmap,
        configs=partcal_config,
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
        config=aoecal_config,
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
        configs=lqcal_config,
        debug_mode=args.debug,
        # gen_plots=bool(args.plot_file),
    )

    msg = f"LQ calibration took {time.time() - start3:.2f} seconds"
    log.info(msg)
    msg = f"Total calibration took {time.time() - start:.2f} seconds"
    log.info(msg)

    if args.plot_file:
        save_dict_to_files(args.plot_file, plot_dicts)

    save_dict_to_files(
        sorted(args.hit_pars),
        {
            tstamp: {
                "pars": {"operations": cal_dicts[tstamp]},
                "results": results_dicts[tstamp],
            }
            for tstamp in cal_dicts
        },
    )

    save_dict_to_files(args.fit_results, object_dicts)
