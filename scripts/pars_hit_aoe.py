import argparse
import json
import logging
import os
import pathlib
import pickle as pkl

import numpy as np
import pygama.pargen.AoE_cal as aoe
from legendmeta import LegendMetadata

argparser = argparse.ArgumentParser()
argparser.add_argument("files", help="files", nargs="*", type=str)
argparser.add_argument("--ecal_file", help="ecal_file", type=str, required=True)
argparser.add_argument("--eres_file", help="eres_file", type=str, required=True)
argparser.add_argument("--inplots", help="in_plot_path", type=str, required=False)

argparser.add_argument("--configs", help="configs", type=str, required=True)
argparser.add_argument("--datatype", help="Datatype", type=str, required=True)
argparser.add_argument("--timestamp", help="Timestamp", type=str, required=True)
argparser.add_argument("--channel", help="Channel", type=str, required=True)

argparser.add_argument("--log", help="log_file", type=str)

argparser.add_argument("--plot_file", help="plot_file", type=str, required=False)
argparser.add_argument("--hit_pars", help="hit_pars", type=str)
argparser.add_argument("--aoe_results", help="aoe_results", type=str)
args = argparser.parse_args()

logging.basicConfig(level=logging.DEBUG, filename=args.log, filemode="w")
logging.getLogger("numba").setLevel(logging.INFO)
logging.getLogger("parse").setLevel(logging.INFO)
logging.getLogger("lgdo").setLevel(logging.INFO)
logging.getLogger("h5py").setLevel(logging.INFO)
logging.getLogger("matplotlib").setLevel(logging.INFO)

configs = LegendMetadata(path=args.configs)
channel_dict = configs.on(args.timestamp, system=args.datatype)["snakemake_rules"][
    "pars_hit_aoecal"
]["inputs"]["aoecal_config"][args.channel]

with open(channel_dict) as r:
    kwarg_dict = json.load(r)

with open(args.ecal_file) as o:
    ecal_dict = json.load(o)
cal_dict = ecal_dict["pars"]
eres_dict = ecal_dict["results"]

with open(args.eres_file, "rb") as o:
    object_dict = pkl.load(o)

if kwarg_dict["run_aoe"] is True:
    kwarg_dict.pop("run_aoe")

    pdf = eval(kwarg_dict.pop("pdf")) if "pdf" in kwarg_dict else aoe.standard_aoe

    sigma_func = (
        eval(kwarg_dict.pop("sigma_func")) if "sigma_func" in kwarg_dict else aoe.sigma_fit
    )

    mean_func = eval(kwarg_dict.pop("mean_func")) if "mean_func" in kwarg_dict else aoe.pol1

    if "plot_options" in kwarg_dict:
        for field, item in kwarg_dict["plot_options"].items():
            kwarg_dict["plot_options"][field]["function"] = eval(item["function"])

    with open(args.files[0]) as f:
        files = f.read().splitlines()
    files = sorted(files)

    try:
        eres = eres_dict[kwarg_dict["cal_energy_param"]]["eres_linear"].copy()

        def eres_func(x):
            return eval(eres["expression"], {"x": x}, eres["pars"])

    except KeyError:

        def eres_func(x):
            return x * np.nan

    cal_dict, out_dict, plot_dict, obj = aoe.aoe_calibration(
        files,
        lh5_path=f"{args.channel}/dsp",
        cal_dicts=cal_dict,
        eres_func=eres_func,
        pdf=pdf,
        mean_func=mean_func,
        sigma_func=sigma_func,
        **kwarg_dict,
    )

    # need to change eres func as can't pickle lambdas
    try:
        obj.eres_func = eres_dict[kwarg_dict["cal_energy_param"]]["eres_linear"].copy()
    except KeyError:
        obj.eres_func = {}
else:
    out_dict = {}
    plot_dict = {}
    obj = None

if args.plot_file:
    common_dict = plot_dict.pop("common") if "common" in list(plot_dict) else None
    if args.inplots:
        with open(args.inplots, "rb") as r:
            out_plot_dict = pkl.load(r)
        out_plot_dict.update(plot_dict)
    else:
        out_plot_dict = plot_dict

    if "common" in list(plot_dict) and common_dict is not None:
        plot_dict("common").update(common_dict)

    pathlib.Path(os.path.dirname(args.plot_file)).mkdir(parents=True, exist_ok=True)
    with open(args.plot_file, "wb") as w:
        pkl.dump(out_plot_dict, w, protocol=pkl.HIGHEST_PROTOCOL)

pathlib.Path(os.path.dirname(args.hit_pars)).mkdir(parents=True, exist_ok=True)
with open(args.hit_pars, "w") as w:
    final_hit_dict = {
        "pars": {"operations": cal_dict},
        "results": {"ecal": eres_dict, "aoe": out_dict},
    }
    json.dump(final_hit_dict, w, indent=4)

pathlib.Path(os.path.dirname(args.aoe_results)).mkdir(parents=True, exist_ok=True)
with open(args.aoe_results, "wb") as w:
    pkl.dump({"ecal": object_dict, "aoe": obj}, w, protocol=pkl.HIGHEST_PROTOCOL)
