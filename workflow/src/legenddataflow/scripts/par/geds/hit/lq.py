from __future__ import annotations

import argparse
import copy
import logging
import pickle as pkl
import time
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from dbetto import TextDB
from dbetto.catalog import Props
from pygama.math.distributions import gaussian
from pygama.pargen.AoE_cal import *  # noqa: F403
from pygama.pargen.lq_cal import *  # noqa: F403
from pygama.pargen.lq_cal import LQCal
from pygama.pargen.utils import load_data

from .....convert_np import convert_dict_np_to_float
from .....log import build_log
from ....pulser_removal import get_pulser_mask

warnings.filterwarnings(action="ignore", category=RuntimeWarning)
log = logging.getLogger(__name__)


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
        debug_mode=debug_mode,
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
    return cal_dicts, get_results_dict(lq), fill_plot_dict(lq, data, plot_options), lq


def run_lq_calibration(
    data,
    cal_dicts,
    results_dicts,
    object_dicts,
    plot_dicts,
    configs,
    debug_mode=False,
    # gen_plots=True,
):
    if isinstance(configs, str | list):
        configs = Props.read_from(configs)

    if configs.pop("run_lq") is True:
        if "plot_options" in configs:
            for field, item in configs["plot_options"].items():
                configs["plot_options"][field]["function"] = eval(item["function"])

        try:
            eres = copy.deepcopy(
                results_dicts[next(iter(results_dicts))]["partition_ecal"][
                    configs["cal_energy_param"]
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
                        configs["cal_energy_param"]
                    ]["eres_linear"]
                )

                def eres_func(x):
                    return eval(eres["expression"], dict(x=x, **eres["parameters"]))

            except KeyError:

                def eres_func(x):
                    return x * np.nan

        log.info("starting lq calibration")
        start = time.time()
        cal_dicts, out_dict, lq_plot_dict, lq_obj = lq_calibration(
            data,
            cal_dicts=cal_dicts,
            energy_param=configs["energy_param"],
            cal_energy_param=configs["cal_energy_param"],
            dt_param=configs["dt_param"],
            eres_func=eres_func,
            cdf=eval(configs.get("cdf", "gaussian")),
            selection_string=f"{configs.pop('cut_field')}&(~is_pulser)",
            plot_options=configs.get("plot_options", None),
            debug_mode=debug_mode | configs.get("debug_mode", False),
        )
        msg = f"lq calibration took {time.time() - start:.2f} seconds"
        log.info(msg)
        # need to change eres func as can't pickle lambdas
        try:
            lq_obj.eres_func = results_dicts[next(iter(results_dicts))][
                "partition_ecal"
            ][configs["cal_energy_param"]]["eres_linear"]
        except KeyError:
            lq_obj.eres_func = {}
    else:
        out_dict = dict.fromkeys(cal_dicts)
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
            plot_dict["common"].update(common_dict)
        elif common_dict is not None:
            plot_dict["common"] = common_dict
        plot_dict.update({"lq": lq_plot_dict})
        out_plot_dicts[tstamp] = plot_dict

    return cal_dicts, out_result_dicts, out_object_dicts, out_plot_dicts


def par_geds_hit_lq() -> None:
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
    argparser.add_argument("--lq-results", help="lq_results", type=str)

    argparser.add_argument("-d", "--debug", help="debug_mode", action="store_true")
    args = argparser.parse_args()

    configs = TextDB(args.configs, lazy=True).on(args.timestamp, system=args.datatype)
    config_dict = configs["snakemake_rules"]["pars_hit_lqcal"]

    build_log(config_dict, args.log)

    channel_dict = config_dict["inputs"]["lqcal_config"][args.channel]
    kwarg_dict = Props.read_from(channel_dict)

    ecal_dict = Props.read_from(args.ecal_file)
    cal_dict = ecal_dict["pars"]["operations"]
    eres_dict = ecal_dict["results"]["ecal"]

    if args.inplots:
        with Path(args.inplots).open("rb") as r:
            plot_dict = pkl.load(r)
    else:
        plot_dict = {}

    with Path(args.eres_file).open("rb") as o:
        object_dict = pkl.load(o)

    if kwarg_dict["run_lq"] is True:
        with Path(args.files[0]).open() as f:
            files = f.read().splitlines()
        files = sorted(files)

        params = [
            "lq80",
            "dt_eff",
            kwarg_dict["energy_param"],
            kwarg_dict["cal_energy_param"],
            kwarg_dict["cut_field"],
        ]

        # load data in
        data, threshold_mask = load_data(
            files,
            args.table_name,
            cal_dict,
            params=params,
            threshold=kwarg_dict.pop("threshold"),
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

        data["run_timestamp"] = args.timestamp

        out_dicts, eres_dicts, plot_dicts, lq_dict = run_lq_calibration(
            data,
            cal_dicts={args.timestamp: cal_dict},
            results_dicts={args.timestamp: eres_dict},
            object_dicts={args.timestamp: object_dict},
            plot_dicts={args.timestamp: plot_dict},
            configs=kwarg_dict,
            debug_mode=args.debug,
        )
        cal_dict = out_dicts[args.timestamp]
        eres_dict = eres_dicts[args.timestamp]
        plot_dict = plot_dicts[args.timestamp]
        lq = lq_dict[args.timestamp]

    else:
        lq = None

    if args.plot_file:
        Path(args.plot_file).parent.mkdir(parents=True, exist_ok=True)
        with Path(args.plot_file).open("wb") as w:
            pkl.dump(plot_dict, w, protocol=pkl.HIGHEST_PROTOCOL)

    final_hit_dict = convert_dict_np_to_float(
        {
            "pars": {"operations": cal_dict},
            "results": dict(**ecal_dict["results"], lq=eres_dict),
        }
    )
    Path(args.hit_pars).parent.mkdir(parents=True, exist_ok=True)
    Props.write_to(args.hit_pars, final_hit_dict)

    final_object_dict = dict(
        **object_dict,
        lq=lq,
    )
    Path(args.lq_results).parent.mkdir(parents=True, exist_ok=True)
    with Path(args.lq_results).open("wb") as w:
        pkl.dump(final_object_dict, w, protocol=pkl.HIGHEST_PROTOCOL)
