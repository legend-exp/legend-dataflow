from __future__ import annotations

import argparse
import copy
import logging
import pickle as pkl
import re
import time
import warnings
from pathlib import Path

import numpy as np
from dbetto import TextDB
from dbetto.catalog import Props
from pygama.pargen.AoE_cal import *  # noqa: F403
from pygama.pargen.AoE_cal import CalAoE
from pygama.pargen.utils import load_data

from .....convert_np import convert_dict_np_to_float
from .....log import build_log
from ....pulser_removal import get_pulser_mask

warnings.filterwarnings(action="ignore", category=RuntimeWarning)

log = logging.getLogger(__name__)


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
    config,
    debug_mode=False,
):
    if isinstance(config, str | list):
        config = Props.read_from(config)

    if config.get("run_aoe", True) is True:
        if "plot_options" in config:
            for field, item in config["plot_options"].items():
                config["plot_options"][field]["function"] = eval(item["function"])

        if "dt_cut" in config and config["dt_cut"] is not None:
            cut_dict = config["dt_cut"]["cut"]
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
                    config["cal_energy_param"]
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
                        config["cal_energy_param"]
                    ]["eres_linear"]
                )

                def eres_func(x):
                    return eval(eres["expression"], dict(x=x, **eres["parameters"]))

            except KeyError:

                def eres_func(x):
                    return x * np.nan

        data["AoE_Uncorr"] = (
            data[config["current_param"]] / data[config["energy_param"]]
        )

        start = time.time()
        log.info("calibrating A/E")

        aoe = CalAoE(
            cal_dicts=cal_dicts,
            cal_energy_param=config["cal_energy_param"],
            eres_func=eres_func,
            pdf=eval(config.get("pdf", "aoe_peak")),
            mean_func=eval(config.get("mean_func", "Pol1")),
            sigma_func=eval(config.get("sigma_func", "SigmaFit")),
            selection_string=f"{config['cut_field']}&(~is_pulser)",
            dt_corr=config.get("dt_corr", False),
            dep_correct=config.get("dep_correct", False),
            dt_cut=config.get("dt_cut", None),
            dt_param=config.get("dt_param", 3),
            high_cut_val=config.get("high_cut_val", 3),
            compt_bands_width=config.get("debug_mode", 20),
            debug_mode=debug_mode | config.get("debug_mode", False),
        )
        aoe.update_cal_dicts(
            {
                "AoE_Uncorr": {
                    "expression": f"{config['current_param']}/{config['energy_param']}",
                    "parameters": {},
                }
            }
        )
        aoe.calibrate(data, "AoE_Uncorr")

        msg = f"A/E calibration completed in {time.time() - start:.2f} seconds"
        log.info(msg)

        out_dict = get_results_dict(aoe)
        aoe_plot_dict = fill_plot_dict(aoe, data, config.get("plot_options", None))

        aoe.pdf = aoe.pdf.name
        # need to change eres func as can't pickle lambdas
        try:
            aoe.eres_func = results_dicts[next(iter(results_dicts))]["partition_ecal"][
                config["cal_energy_param"]
            ]["eres_linear"]
        except KeyError:
            aoe.eres_func = {}
    else:
        out_dict = dict.fromkeys(cal_dicts)
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


def par_geds_hit_aoe() -> None:
    argparser = argparse.ArgumentParser()
    argparser.add_argument("files", help="files", nargs="*", type=str)
    argparser.add_argument(
        "--pulser-file", help="pulser_file", type=str, required=False
    )
    argparser.add_argument(
        "--tcm-filelist", help="tcm_filelist", type=str, required=False
    )

    argparser.add_argument("--ecal-file", help="ecal_file", type=str, required=True)
    argparser.add_argument("--eres-file", help="eres_file", type=str, required=True)
    argparser.add_argument("--inplots", help="in_plot_path", type=str, required=False)

    argparser.add_argument("--configs", help="configs", type=str, required=True)
    argparser.add_argument("--log", help="log_file", type=str)

    argparser.add_argument("--datatype", help="Datatype", type=str, required=True)
    argparser.add_argument("--timestamp", help="Timestamp", type=str, required=True)
    argparser.add_argument("--channel", help="Channel", type=str, required=True)
    argparser.add_argument("--table-name", help="table name", type=str, required=True)

    argparser.add_argument("--plot-file", help="plot_file", type=str, required=False)
    argparser.add_argument("--hit-pars", help="hit_pars", type=str)
    argparser.add_argument("--aoe-results", help="aoe_results", type=str)

    argparser.add_argument("-d", "--debug", help="debug_mode", action="store_true")
    args = argparser.parse_args()

    configs = TextDB(args.configs, lazy=True).on(args.timestamp, system=args.datatype)
    config_dict = configs["snakemake_rules"]["pars_hit_aoecal"]

    build_log(config_dict, args.log)

    channel_dict = config_dict["inputs"]["aoecal_config"][args.channel]
    kwarg_dict = Props.read_from(channel_dict)

    ecal_dict = Props.read_from(args.ecal_file)
    cal_dict = ecal_dict["pars"]
    eres_dict = ecal_dict["results"]["ecal"]

    with Path(args.eres_file).open("rb") as o:
        object_dict = pkl.load(o)

    if args.inplots:
        with Path(args.inplots).open("rb") as r:
            out_plot_dict = pkl.load(r)
    else:
        out_plot_dict = {}

    with Path(args.files[0]).open() as f:
        files = sorted(f.read().splitlines())

    if kwarg_dict["run_aoe"] is True:
        params = [
            kwarg_dict["current_param"],
            "tp_0_est",
            "tp_99",
            kwarg_dict["energy_param"],
            kwarg_dict["cal_energy_param"],
            kwarg_dict["cut_field"],
            "timestamp",
        ]

        if "dt_param" in kwarg_dict:
            params += kwarg_dict["dt_param"]
        else:
            params.append("dt_eff")

        if "dt_cut" in kwarg_dict and kwarg_dict["dt_cut"] is not None:
            cal_dict.update(kwarg_dict["dt_cut"]["cut"])
            params.append(kwarg_dict["dt_cut"]["out_param"])

        # load data in
        data, threshold_mask = load_data(
            files,
            args.table_name,
            cal_dict,
            params=params,
            threshold=kwarg_dict["threshold"],
            return_selection_mask=True,
        )
        msg = f"Loaded {len(data)} events"
        log.info(msg)

        mask = get_pulser_mask(
            pulser_file=args.pulser_file,
        )

        data["is_pulser"] = mask[threshold_mask]

        msg = f"{len(data.query('~is_pulser'))}  non pulser events"
        log.info(msg)

        with Path(args.files[0]).open() as f:
            files = f.read().splitlines()
        files = sorted(files)

        data["run_timestamp"] = args.timestamp

        cal_dict, results_dicts, object_dicts, plot_dicts = run_aoe_calibration(
            data,
            {args.timestamp: cal_dict},
            {args.timestamp: eres_dict},
            {args.timestamp: object_dict},
            {args.timestamp: out_plot_dict},
            kwarg_dict,
            debug_mode=args.debug,
        )
        cal_dict = cal_dict[args.timestamp]
        results_dict = results_dicts[args.timestamp]
        aoe = object_dicts[args.timestamp]
        plot_dict = plot_dicts[args.timestamp]
    else:
        aoe = None
        plot_dict = out_plot_dict

    if args.plot_file:
        common_dict = plot_dict.pop("common") if "common" in list(plot_dict) else None
        out_plot_dict.update({"aoe": plot_dict})

        if "common" in list(out_plot_dict) and common_dict is not None:
            out_plot_dict["common"].update(common_dict)
        elif common_dict is not None:
            out_plot_dict["common"] = common_dict

        Path(args.plot_file).parent.mkdir(parents=True, exist_ok=True)
        with Path(args.plot_file).open("wb") as w:
            pkl.dump(out_plot_dict, w, protocol=pkl.HIGHEST_PROTOCOL)

    Path(args.hit_pars).parent.mkdir(parents=True, exist_ok=True)
    final_hit_dict = {
        "pars": {"operations": cal_dict},
        "results": dict(**ecal_dict["results"], aoe=results_dict),
    }

    final_hit_dict = convert_dict_np_to_float(final_hit_dict)

    Props.write_to(args.hit_pars, final_hit_dict)

    Path(args.aoe_results).parent.mkdir(parents=True, exist_ok=True)
    with Path(args.aoe_results).open("wb") as w:
        pkl.dump(dict(**object_dict, aoe=aoe), w, protocol=pkl.HIGHEST_PROTOCOL)
