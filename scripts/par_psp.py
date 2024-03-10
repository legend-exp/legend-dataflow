import argparse
import json
import os
import pathlib
from legendmeta.catalog import Props
from util.FileKey import ChannelProcKey
import numpy as np
import matplotlib.pyplot as pyplot
import matplotlib as mpl 
from datetime import datetime
import pickle as pkl
mpl.use("Agg")


argparser = argparse.ArgumentParser()
argparser.add_argument("--input", help="input files", nargs="*", type=str, required=True)
argparser.add_argument("--output", help="output file", nargs="*", type=str, required=True)
argparser.add_argument("--in_plots", help="input plot files", nargs="*", type=str, required=False)
argparser.add_argument("--out_plots", help="output plot files", nargs="*", type=str, required=False)
argparser.add_argument("--in_obj", help="input object files", nargs="*", type=str, required=False)
argparser.add_argument("--out_obj", help="output object files", nargs="*", type=str, required=False)
args = argparser.parse_args()

conf = LegendMetadata(path=args.configs)
configs = conf.on(args.timestamp, system=args.datatype)
merge_config = configs["snakemake_rules"]["pars_psp"]["inputs"]["config"][
    args.channel
]

ave_fields = merge_config["average_fields"]

# partitions could be different for different channels - do separately for each channel
in_dicts = {}
for file in args.input:
    tstamp = ChannelProcKey.get_filekey_from_pattern(os.path.basename(file)).timestamp
    in_dicts[tstamp] = Props.read_from(file)

plot_dict = {}
for field in ave_fields:
    keys = field.split(".")
    vals = []
    for tstamp in in_dicts:
        val = in_dicts[tstamp]
        for key in keys:
            val = val[key]
        vals.append(val)
    if len(vals[~np.isnan(vals)]) ==0:
        mean = np.nan 
    else:
        mean = np.nanmean(vals)
    for tstamp in in_dicts:
        val = in_dicts[tstamp]
        for key in keys:
            val = val[key]
        val = mean

    fig = plt.figure()
    plt.scatter([datetime.strptime(tstamp,'%Y%m%dT%H%M%SZ') for tstamp in in_dicts] , vals)
    plt.axhline(y=mean, color='r', linestyle='-')
    plt.xlabel("time")
    plt.ylabel("value")
    plt.title(f"{field} over time")
    plot_dict[field] = fig
    plt.close()

for file in args.output:
    tstamp = ChannelProcKey.get_filekey_from_pattern(os.path.basename(file)).timestamp
    with open(file, "w") as f:
        json.dump(in_dicts[tstamp], f, indent=2)


if args.out_plots:
    for file in args.out_plots:
        tstamp = ChannelProcKey.get_filekey_from_pattern(os.path.basename(file)).timestamp
        if args.in_plots:
            for infile in args.in_plots:
                if tstamp in infile:
                    with open(infile, "rb") as f:
                        old_plot_dict = pkl.load(f)
                    break
            new_plot_dict = old_plot_dict.update({"psp": plot_dict})
        else:
            new_plot_dict = {"psp": plot_dict}
        with open(file, "w") as f:
            pkl.dump(new_plot_dict, file, protocol=pkl.HIGHEST_PROTOCOL)

if args.out_obj:
    for file in args.out_obj:
        tstamp = ChannelProcKey.get_filekey_from_pattern(os.path.basename(file)).timestamp
        if args.in_obj:
            for infile in args.in_obj:
                if tstamp in infile:
                    with open(infile, "rb") as f:
                        old_obj_dict = pkl.load(f)
                    break
            new_obj_dict = old_obj_dict
        else:
            new_obj_dict = {}
        with open(file, "w") as f:
            pkl.dump(new_obj_dict, file, protocol=pkl.HIGHEST_PROTOCOL)