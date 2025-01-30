import argparse
import pickle as pkl
from datetime import datetime
from pathlib import Path

import matplotlib as mpl
import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import numpy as np
from dbetto.catalog import Props
from legendmeta import LegendMetadata

from ..FileKey import ChannelProcKey

mpl.use("Agg")


argparser = argparse.ArgumentParser()
argparser.add_argument(
    "--input", help="input files", nargs="*", type=str, required=True
)
argparser.add_argument(
    "--output", help="output file", nargs="*", type=str, required=True
)
argparser.add_argument(
    "--in_plots", help="input plot files", nargs="*", type=str, required=False
)
argparser.add_argument(
    "--out_plots", help="output plot files", nargs="*", type=str, required=False
)
argparser.add_argument(
    "--in_obj", help="input object files", nargs="*", type=str, required=False
)
argparser.add_argument(
    "--out_obj", help="output object files", nargs="*", type=str, required=False
)

argparser.add_argument("--log", help="log_file", type=str)
argparser.add_argument("--configs", help="configs", type=str, required=True)

argparser.add_argument("--datatype", help="Datatype", type=str, required=True)
argparser.add_argument("--timestamp", help="Timestamp", type=str, required=True)
argparser.add_argument("--channel", help="Channel", type=str, required=True)
args = argparser.parse_args()

configs = LegendMetadata(args.configs, lazy=True).on(
    args.timestamp, system=args.datatype
)
merge_config = Props.read_from(
    configs["snakemake_rules"]["pars_psp"]["inputs"]["psp_config"][args.channel]
)

ave_fields = merge_config["average_fields"]

# partitions could be different for different channels - do separately for each channel
in_dicts = {}
for file in args.input:
    tstamp = ChannelProcKey.get_filekey_from_pattern(Path(file).name).timestamp
    in_dicts[tstamp] = Props.read_from(file)

plot_dict = {}
for field in ave_fields:
    keys = field.split(".")
    vals = []
    for _tstamp, tstamp_dict in in_dicts.items():
        val = tstamp_dict.copy()
        for key in keys:
            val = val[key]
        vals.append(val)
        if "dsp" in tstamp_dict:
            tmp_dict = tstamp_dict["dsp"]
        else:
            tmp_dict = {}
            tstamp_dict["dsp"] = tmp_dict
        for i, key in enumerate(keys):
            if i == len(keys) - 1:
                tmp_dict[key] = val
            else:
                if key in tmp_dict:
                    tmp_dict = tmp_dict[key]
                else:
                    tmp_dict[key] = {}
                    tmp_dict = tmp_dict[key]
    if isinstance(vals[0], str):
        if "*" in vals[0]:
            unit = vals[0].split("*")[1]
            rounding = len(val.split("*")[0].split(".")[-1]) if "." in vals[0] else 16
            vals = np.array([float(val.split("*")[0]) for val in vals])
        else:
            unit = None
            rounding = 16
    else:
        vals = np.array(vals)
        unit = None
        rounding = 16

    mean_val = np.nan if len(vals[~np.isnan(vals)]) == 0 else np.nanmedian(vals)
    mean = f"{round(mean_val, rounding)}*{unit}" if unit is not None else mean_val

    for _tstamp, tstamp_dict in in_dicts.items():
        val = tstamp_dict
        for i, key in enumerate(keys):
            if i == len(keys) - 1:
                val[key] = mean
            else:
                val = val[key]

    fig = plt.figure()
    plt.scatter(
        [datetime.strptime(tstamp, "%Y%m%dT%H%M%SZ") for tstamp in in_dicts], vals
    )
    plt.axhline(y=mean_val, color="r", linestyle="-")
    plt.xlabel("time")
    if unit is not None:
        plt.ylabel(f"value {unit}")
    else:
        plt.ylabel("value")
    plt.gca().xaxis.set_major_formatter(mdates.DateFormatter("%d/%m/%y"))
    plt.gcf().autofmt_xdate()
    plt.title(field)
    plot_dict[field] = fig
    plt.close()

for file in args.output:
    tstamp = ChannelProcKey.get_filekey_from_pattern(Path(file).name).timestamp
    Props.write_to(file, in_dicts[tstamp])

if args.out_plots:
    for file in args.out_plots:
        tstamp = ChannelProcKey.get_filekey_from_pattern(Path(file).name).timestamp
        if args.in_plots:
            for infile in args.in_plots:
                if tstamp in infile:
                    with Path(infile).open("rb") as f:
                        old_plot_dict = pkl.load(f)
                    break
            old_plot_dict.update({"psp": plot_dict})
            new_plot_dict = old_plot_dict
        else:
            new_plot_dict = {"psp": plot_dict}
        with Path(file).open("wb") as f:
            pkl.dump(new_plot_dict, f, protocol=pkl.HIGHEST_PROTOCOL)

if args.out_obj:
    for file in args.out_obj:
        tstamp = ChannelProcKey.get_filekey_from_pattern(Path(file).name).timestamp
        if args.in_obj:
            for infile in args.in_obj:
                if tstamp in infile:
                    with Path(infile).open("rb") as f:
                        old_obj_dict = pkl.load(f)
                    break
            new_obj_dict = old_obj_dict
        else:
            new_obj_dict = {}
        with Path(file).open("wb") as f:
            pkl.dump(new_obj_dict, f, protocol=pkl.HIGHEST_PROTOCOL)
