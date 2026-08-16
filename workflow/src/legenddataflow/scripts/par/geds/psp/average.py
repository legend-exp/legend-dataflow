"""Console script ``par-geds-psp-average``: average run-level DSP parameters
over a partition to produce psp parameters, with stability plots."""

from __future__ import annotations

import argparse
import pickle as pkl
from datetime import datetime
from pathlib import Path

import matplotlib as mpl
import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import numpy as np
from dbetto import AttrsDict, Props
from legenddataflowscripts.utils import (
    check_input_files,
    convert_dict_np_to_float,
    get_channel_config,
    get_rule_config,
    prepare_output_paths,
)

from legenddataflow.methods import ChannelProcKey

mpl.use("Agg")


def par_geds_psp_average() -> None:
    argparser = argparse.ArgumentParser()
    argparser.add_argument(
        "--input", help="input files", nargs="*", type=str, required=True
    )
    argparser.add_argument(
        "--output", help="output file", nargs="*", type=str, required=True
    )
    argparser.add_argument(
        "--in-plots", help="input plot files", nargs="*", type=str, required=False
    )
    argparser.add_argument(
        "--out-plots", help="output plot files", nargs="*", type=str, required=False
    )
    argparser.add_argument(
        "--in-obj", help="input object files", nargs="*", type=str, required=False
    )
    argparser.add_argument(
        "--out-obj", help="output object files", nargs="*", type=str, required=False
    )

    argparser.add_argument("--log", help="log_file", type=str)
    argparser.add_argument("--configs", help="configs", type=str, required=True)

    argparser.add_argument("--datatype", help="Datatype", type=str, required=True)
    argparser.add_argument("--timestamp", help="Timestamp", type=str, required=True)
    argparser.add_argument("--channel", help="Channel", type=str, required=True)
    args = argparser.parse_args()

    config_dict = get_rule_config(
        args.configs, "pars_psp", args.timestamp, args.datatype
    )
    merge_config = Props.read_from(
        get_channel_config(
            config_dict["inputs"]["psp_config"], args.channel, name="psp_config"
        )
    )

    ave_fields = merge_config["average_fields"]

    check_input_files(args.input, "--input")
    prepare_output_paths(
        *(args.output or []), *(args.out_plots or []), *(args.out_obj or [])
    )

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
                elif key in tmp_dict:
                    tmp_dict = tmp_dict[key]
                else:
                    tmp_dict[key] = {}
                    tmp_dict = tmp_dict[key]
        if isinstance(vals[0], str):
            try:
                if "*" in vals[0]:
                    unit = vals[0].split("*")[1]
                    rounding = (
                        len(vals[0].split("*")[0].split(".")[-1])
                        if "." in vals[0]
                        else 16
                    )
                    vals = np.array([float(val.split("*")[0]) for val in vals])
                else:
                    unit = None
                    rounding = 16
                    vals = np.array([float(val) for val in vals])
            except (ValueError, IndexError) as err:
                msg = (
                    f"cannot average field {field!r}: values {vals} are not "
                    "numbers or 'number*unit' strings"
                )
                raise ValueError(msg) from err
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
        Props.write_to(
            file, AttrsDict(convert_dict_np_to_float(in_dicts[tstamp])).to_dict()
        )

    def _find_input_for_timestamp(infiles, tstamp, outfile):
        for infile in infiles:
            if tstamp in infile:
                with Path(infile).open("rb") as f:
                    return pkl.load(f)
        msg = (
            f"no input file matching timestamp {tstamp} (from output {outfile}) "
            f"found among {infiles}"
        )
        raise ValueError(msg)

    if args.out_plots:
        for file in args.out_plots:
            tstamp = ChannelProcKey.get_filekey_from_pattern(Path(file).name).timestamp
            if args.in_plots:
                new_plot_dict = _find_input_for_timestamp(args.in_plots, tstamp, file)
                new_plot_dict.update({"psp": plot_dict})
            else:
                new_plot_dict = {"psp": plot_dict}
            with Path(file).open("wb") as f:
                pkl.dump(new_plot_dict, f, protocol=pkl.HIGHEST_PROTOCOL)

    if args.out_obj:
        for file in args.out_obj:
            tstamp = ChannelProcKey.get_filekey_from_pattern(Path(file).name).timestamp
            if args.in_obj:
                new_obj_dict = _find_input_for_timestamp(args.in_obj, tstamp, file)
            else:
                new_obj_dict = {}
            with Path(file).open("wb") as f:
                pkl.dump(new_obj_dict, f, protocol=pkl.HIGHEST_PROTOCOL)
