import argparse
import json
import logging
import os
import pathlib
import pickle as pkl

import numpy as np
from legendmeta import LegendMetadata
from pygama.pargen.AoE_cal import *
from pygama.pargen.ecal_th import *
from util.FileKey import *

argparser = argparse.ArgumentParser()
argparser.add_argument("--input_files", help="files", type=str, nargs="*", required=True)
argparser.add_argument("--ecal_file", help="ecal_file", type=str, nargs="*", required=True)
argparser.add_argument("--eres_file", help="eres_file", type=str, nargs="*", required=True)
argparser.add_argument("--inplots", help="eres_file", type=str, nargs="*", required=True)

argparser.add_argument("--configs", help="configs", type=str, required=True)
argparser.add_argument("--timestamp", help="Datatype", type=str, required=True)
argparser.add_argument("--datatype", help="Datatype", type=str, required=True)
argparser.add_argument("--channel", help="Channel", type=str, required=True)

argparser.add_argument("--log", help="log_file", type=str)

argparser.add_argument("--plot_file", help="plot_file", type=str, nargs="*", required=False)
argparser.add_argument("--hit_pars", help="hit_pars", nargs="*", type=str)
argparser.add_argument("--aoe_results", help="aoe_results", nargs="*", type=str)
args = argparser.parse_args()

logging.basicConfig(level=logging.DEBUG, filename=args.log, filemode="w")
logging.getLogger("numba").setLevel(logging.INFO)
logging.getLogger("parse").setLevel(logging.INFO)
logging.getLogger("lgdo").setLevel(logging.INFO)
logging.getLogger("h5py").setLevel(logging.INFO)
logging.getLogger("matplotlib").setLevel(logging.INFO)


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


configs = LegendMetadata(path=args.configs)
channel_dict = configs.on(args.timestamp, system=args.datatype)["snakemake_rules"]["pars_pht"][
    "inputs"
]["par_pht_config"][args.channel]

with open(channel_dict) as r:
    kwarg_dict = json.load(r)

cal_dict = {}
results_dicts = {}
if isinstance(args.ecal_file, list):
    for ecal in args.ecal_file:
        with open(ecal) as o:
            cal = json.load(o)

        fk = ChannelProcKey.get_filekey_from_pattern(os.path.basename(ecal))
        cal_dict[fk.timestamp] = cal["pars"]
        results_dicts[fk.timestamp] = cal["results"]
else:
    with open(args.ecal_file) as o:
        cal = json.load(o)

    fk = ChannelProcKey.get_filekey_from_pattern(os.path.basename(args.ecal_file))
    cal_dict[fk.timestamp] = cal["pars"]
    results_dicts[fk.timestamp] = cal["results"]

object_dict = {}
if isinstance(args.eres_file, list):
    for ecal in args.eres_file:
        with open(ecal, "rb") as o:
            cal = pkl.load(o)
        fk = ChannelProcKey.get_filekey_from_pattern(os.path.basename(ecal))
        object_dict[fk.timestamp] = cal
else:
    with open(args.eres_file, "rb") as o:
        cal = pkl.load(o)
    fk = ChannelProcKey.get_filekey_from_pattern(os.path.basename(args.eres_file))
    object_dict[fk.timestamp] = cal

inplots_dict = {}
if args.inplots:
    if isinstance(args.inplots, list):
        for ecal in args.inplots:
            with open(ecal, "rb") as o:
                cal = pkl.load(o)
            fk = ChannelProcKey.get_filekey_from_pattern(os.path.basename(ecal))
            inplots_dict[fk.timestamp] = cal
    else:
        with open(args.inplots, "rb") as o:
            cal = pkl.load(o)
        fk = ChannelProcKey.get_filekey_from_pattern(os.path.basename(args.inplots))
        inplots_dict[fk.timestamp] = cal

ecal_options = kwarg_dict.pop("ecal")
aoe_options = kwarg_dict.pop("aoe")


if "plot_options" in ecal_options:
    for field, item in ecal_options["plot_options"].items():
        ecal_options["plot_options"][field]["function"] = eval(item["function"])

if "plot_options" in aoe_options:
    for field, item in aoe_options["plot_options"].items():
        aoe_options["plot_options"][field]["function"] = eval(item["function"])


# sort files in dictionary where keys are first timestamp from run
if isinstance(args.input_files, list):
    files = []
    for file in args.input_files:
        with open(file) as f:
            files += f.read().splitlines()
else:
    with open(args.input_files) as f:
        files = f.read().splitlines()

files = sorted(
    np.unique(files)
)  # need this as sometimes files get double counted as it somehow puts in the p%-* filelist and individual runs also

final_dict = {}
all_file = run_splitter(sorted(files))
for filelist in all_file:
    fk = ProcessingFileKey.get_filekey_from_pattern(os.path.basename(sorted(filelist)[0]))
    timestamp = fk.timestamp
    final_dict[timestamp] = sorted(filelist)

# run energy supercal
ecal_results, ecal_plots, ecal_obj = partition_energy_cal_th(
    final_dict, lh5_path=f"{args.channel}/dsp", hit_dict=cal_dict, **ecal_options
)

# run aoe cal
if aoe_options.pop("run_aoe") is True:
    pdf = eval(aoe_options.pop("pdf")) if "pdf" in aoe_options else standard_aoe

    if "mean_func" in aoe_options:
        mean_func = eval(aoe_options.pop("mean_func"))
    else:
        mean_func = pol1

    if "sigma_func" in aoe_options:
        sigma_func = eval(aoe_options.pop("sigma_func"))
    else:
        sigma_func = sigma_fit

    try:
        eres = ecal_results[aoe_options["cal_energy_param"]]["eres_linear"].copy()

        def eres_func(x):
            return eval(eres["expression"], {"x": x}, eres["parameters"])

        if np.isnan(eres_func(2000)):
            raise RuntimeError
    except:
        try:
            eres = results_dicts[list(results_dicts)[0]][aoe_options["cal_energy_param"]][
                "eres_linear"
            ].copy()

            def eres_func(x):
                return eval(eres["expression"], {"x": x}, eres["parameters"])

        except:

            def eres_func():
                return np.nan

    cal_dict, out_dict, plot_dict, aoe_obj = aoe_calibration(
        final_dict,
        lh5_path=f"{args.channel}/dsp",
        cal_dicts=cal_dict,
        eres_func=eres_func,
        pdf=pdf,
        mean_func=mean_func,
        sigma_func=sigma_func,
        **aoe_options,
    )

    # need to change eres func as can't pickle lambdas
    try:
        aoe_obj.eres_func = eres_dict[kwarg_dict["cal_energy_param"]]["eres_linear"].copy()
    except:
        aoe_obj.eres_func = {}
else:
    out_dict = {}
    plot_dict = {}
    aoe_obj = None


if args.plot_file:
    common_dict = plot_dict.pop("common") if "common" in list(plot_dict) else None

    if isinstance(args.plot_file, list):
        for plot_file in args.plot_file:
            fk = ChannelProcKey.get_filekey_from_pattern(os.path.basename(plot_file))
            if args.inplots:
                out_plot_dict = inplots_dict[fk.timestamp]
                out_plot_dict.update({"aoe": plot_dict})
                out_plot_dict.update({"partition_ecal": ecal_plots})
            else:
                out_plot_dict = {"aoe": plot_dict}
                out_plot_dict.update({"partition_ecal": ecal_plots})

            if "common" in list(plot_dict) and common_dict is not None:
                plot_dict("common").update(common_dict)

            pathlib.Path(os.path.dirname(plot_file)).mkdir(parents=True, exist_ok=True)
            with open(plot_file, "wb") as w:
                pkl.dump(out_plot_dict, w, protocol=pkl.HIGHEST_PROTOCOL)
    else:
        if args.inplots:
            fk = ChannelProcKey.get_filekey_from_pattern(os.path.basename(args.plot_file))
            out_plot_dict = inplots_dict[fk.timestamp]
            out_plot_dict.update(plot_dict)
            out_plot_dict.update({"partition_ecal": ecal_plots})
        else:
            out_plot_dict = plot_dict
        pathlib.Path(os.path.dirname(args.plot_file)).mkdir(parents=True, exist_ok=True)
        with open(out, "wb") as w:
            pkl.dump(out_plot_dict, w, protocol=pkl.HIGHEST_PROTOCOL)


for out in sorted(args.hit_pars):
    fk = ChannelProcKey.get_filekey_from_pattern(os.path.basename(out))
    final_hit_dict = {
        "pars": {"operations": cal_dict[fk.timestamp]},
        "results": {
            "ecal": results_dicts[fk.timestamp],
            "partition_ecal": ecal_results,
            "aoe": out_dict,
        },
    }
    pathlib.Path(os.path.dirname(out)).mkdir(parents=True, exist_ok=True)
    with open(out, "w") as w:
        json.dump(final_hit_dict, w, indent=4)

for out in args.aoe_results:
    fk = ChannelProcKey.get_filekey_from_pattern(os.path.basename(out))
    final_object_dict = {
        "ecal": object_dict[fk.timestamp],
        "partition_ecal": ecal_obj,
        "aoe": aoe_obj,
    }
    pathlib.Path(os.path.dirname(out)).mkdir(parents=True, exist_ok=True)
    with open(out, "wb") as w:
        pkl.dump(final_object_dict, w, protocol=pkl.HIGHEST_PROTOCOL)
