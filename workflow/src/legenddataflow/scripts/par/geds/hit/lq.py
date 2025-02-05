from __future__ import annotations

import argparse
import pickle as pkl
import warnings
from pathlib import Path

import numpy as np
from dbetto import TextDB
from dbetto.catalog import Props
from legendmeta import LegendMetadata
from pygama.math.distributions import gaussian
from pygama.pargen.AoE_cal import *  # noqa: F403
from pygama.pargen.data_cleaning import get_tcm_pulser_ids
from pygama.pargen.lq_cal import *  # noqa: F403
from pygama.pargen.lq_cal import LQCal
from pygama.pargen.utils import load_data

from ....convert_np import convert_dict_np_to_float
from ....log import build_log

warnings.filterwarnings(action="ignore", category=RuntimeWarning)


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
    if plot_dict is not None:
        for key, item in plot_options.items():
            if item["options"] is not None:
                plot_dict[key] = item["function"](lq_class, data, **item["options"])
            else:
                plot_dict[key] = item["function"](lq_class, data)
    else:
        plot_dict = {}
    return plot_dict


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
    argparser.add_argument("--metadata", help="metadata", type=str, required=True)
    argparser.add_argument("--log", help="log_file", type=str)

    argparser.add_argument("--datatype", help="Datatype", type=str, required=True)
    argparser.add_argument("--timestamp", help="Timestamp", type=str, required=True)
    argparser.add_argument("--channel", help="Channel", type=str, required=True)

    argparser.add_argument("--plot-file", help="plot_file", type=str, required=False)
    argparser.add_argument("--hit-pars", help="hit_pars", type=str)
    argparser.add_argument("--lq-results", help="lq_results", type=str)

    argparser.add_argument("-d", "--debug", help="debug_mode", action="store_true")
    args = argparser.parse_args()

    configs = TextDB(args.configs, lazy=True).on(args.timestamp, system=args.datatype)
    config_dict = configs["snakemake_rules"]["pars_hit_lqcal"]

    log = build_log(config_dict, args.log)

    meta = LegendMetadata(path=args.metadata)
    channel_dict = meta.channelmap(args.timestamp, system=args.datatype)
    channel = f"ch{channel_dict[args.channel].daq.rawid:07}"

    channel_dict = config_dict["inputs"]["lqcal_config"][args.channel]
    kwarg_dict = Props.read_from(channel_dict)

    ecal_dict = Props.read_from(args.ecal_file)
    cal_dict = ecal_dict["pars"]["operations"]
    eres_dict = ecal_dict["results"]["ecal"]

    with Path(args.eres_file).open("rb") as o:
        object_dict = pkl.load(o)

    if kwarg_dict["run_lq"] is True:
        kwarg_dict.pop("run_lq")

        cdf = eval(kwarg_dict.pop("cdf")) if "cdf" in kwarg_dict else gaussian

        if "plot_options" in kwarg_dict:
            for field, item in kwarg_dict["plot_options"].items():
                kwarg_dict["plot_options"][field]["function"] = eval(item["function"])

        with Path(args.files[0]).open() as f:
            files = f.read().splitlines()
        files = sorted(files)

        try:
            eres = eres_dict[kwarg_dict["cal_energy_param"]]["eres_linear"].copy()

            def eres_func(x):
                return eval(eres["expression"], dict(x=x, **eres["parameters"]))

        except KeyError:

            def eres_func(x):
                return x * np.nan

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
            f"{channel}/dsp",
            cal_dict,
            params=params,
            threshold=kwarg_dict.pop("threshold"),
            return_selection_mask=True,
        )

        if args.pulser_file:
            pulser_dict = Props.read_from(args.pulser_file)
            mask = np.array(pulser_dict["mask"])
            if "pulser_multiplicity_threshold" in kwarg_dict:
                kwarg_dict.pop("pulser_multiplicity_threshold")

        elif args.tcm_filelist:
            # get pulser mask from tcm files
            with Path(args.tcm_filelist).open() as f:
                tcm_files = f.read().splitlines()
            tcm_files = sorted(np.unique(tcm_files))
            ids, mask = get_tcm_pulser_ids(
                tcm_files, channel, kwarg_dict.pop("pulser_multiplicity_threshold")
            )
        else:
            msg = "No pulser file or tcm filelist provided"
            raise ValueError(msg)

        data["is_pulser"] = mask[threshold_mask]

        lq = LQCal(
            cal_dict,
            kwarg_dict["cal_energy_param"],
            kwarg_dict["dt_param"],
            eres_func,
            cdf,
            selection_string=f"{kwarg_dict.pop('cut_field')}&(~is_pulser)",
            debug_mode=args.debug_mode | kwarg_dict.get("debug_mode", False),
        )

        data["LQ_Ecorr"] = np.divide(data["lq80"], data[kwarg_dict["energy_param"]])

        lq.update_cal_dicts(
            {
                "LQ_Ecorr": {
                    "expression": f"lq80/{kwarg_dict['energy_param']}",
                    "parameters": {},
                }
            }
        )

        lq.calibrate(data, "LQ_Ecorr")
        log.info("Calibrated LQ")

        out_dict = get_results_dict(lq)
        plot_dict = fill_plot_dict(lq, data, kwarg_dict.get("plot_options", None))

        # need to change eres func as can't pickle lambdas
        try:
            lq.eres_func = eres_dict[kwarg_dict["cal_energy_param"]][
                "eres_linear"
            ].copy()
        except KeyError:
            lq.eres_func = {}
    else:
        out_dict = {}
        plot_dict = {}
        lq = None

    if args.plot_file:
        common_dict = plot_dict.pop("common") if "common" in list(plot_dict) else None
        if args.inplots:
            with Path(args.inplots).open("rb") as r:
                out_plot_dict = pkl.load(r)
            out_plot_dict.update({"lq": plot_dict})
        else:
            out_plot_dict = {"lq": plot_dict}

        if "common" in list(out_plot_dict) and common_dict is not None:
            out_plot_dict["common"].update(common_dict)
        elif common_dict is not None:
            out_plot_dict["common"] = common_dict

        Path(args.plot_file).parent.mkdir(parents=True, exist_ok=True)
        with Path(args.plot_file).open("wb") as w:
            pkl.dump(out_plot_dict, w, protocol=pkl.HIGHEST_PROTOCOL)

    final_hit_dict = convert_dict_np_to_float(
        {
            "pars": {"operations": cal_dict},
            "results": dict(**eres_dict, lq=out_dict),
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
