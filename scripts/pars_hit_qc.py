from __future__ import annotations

import argparse
import json
import logging
import os
import pathlib
import pickle as pkl
import re
import warnings

os.environ["PYGAMA_PARALLEL"] = "false"
os.environ["PYGAMA_FASTMATH"] = "false"

import numpy as np
from legendmeta import LegendMetadata
from legendmeta.catalog import Props
from lgdo.lh5 import ls
from pygama.pargen.data_cleaning import (
    generate_cut_classifiers,
    get_keys,
    get_tcm_pulser_ids,
)
from pygama.pargen.utils import load_data

log = logging.getLogger(__name__)

warnings.filterwarnings(action="ignore", category=RuntimeWarning)


if __name__ == "__main__":
    argparser = argparse.ArgumentParser()
    argparser.add_argument("--cal_files", help="cal_files", nargs="*", type=str)
    argparser.add_argument("--fft_files", help="fft_files", nargs="*", type=str)
    argparser.add_argument("--tcm_filelist", help="tcm_filelist", type=str, required=False)
    argparser.add_argument("--pulser_file", help="pulser_file", type=str, required=False)

    argparser.add_argument("--configs", help="config", type=str, required=True)
    argparser.add_argument("--datatype", help="Datatype", type=str, required=True)
    argparser.add_argument("--timestamp", help="Timestamp", type=str, required=True)
    argparser.add_argument("--channel", help="Channel", type=str, required=True)
    argparser.add_argument("--tier", help="tier", type=str, default="hit")

    argparser.add_argument("--log", help="log_file", type=str)

    argparser.add_argument("--plot_path", help="plot_path", type=str, required=False)
    argparser.add_argument("--save_path", help="save_path", type=str)
    args = argparser.parse_args()

    logging.basicConfig(level=logging.DEBUG, filename=args.log, filemode="w")
    logging.getLogger("numba").setLevel(logging.INFO)
    logging.getLogger("parse").setLevel(logging.INFO)
    logging.getLogger("lgdo").setLevel(logging.INFO)
    logging.getLogger("h5py").setLevel(logging.INFO)
    logging.getLogger("matplotlib").setLevel(logging.INFO)
    logging.getLogger("legendmeta").setLevel(logging.INFO)

    # get metadata dictionary
    configs = LegendMetadata(path=args.configs)
    channel_dict = configs.on(args.timestamp, system=args.datatype)["snakemake_rules"]
    channel_dict = channel_dict["pars_hit_qc"]["inputs"]["qc_config"][args.channel]

    kwarg_dict = Props.read_from(channel_dict)

    kwarg_dict_fft = kwarg_dict["fft_fields"]
    if len(args.fft_files) > 0:
        fft_fields = get_keys(
            [
                key.replace(f"{args.channel}/dsp/", "")
                for key in ls(args.fft_files[0], f"{args.channel}/dsp/")
            ],
            kwarg_dict_fft["cut_parameters"],
        )

        fft_data = load_data(
            args.fft_files,
            f"{args.channel}/dsp",
            {},
            [*fft_fields, "timestamp", "trapTmax"],
        )

        discharges = fft_data["t_sat_lo"] > 0
        discharge_timestamps = np.where(fft_data["timestamp"][discharges])[0]
        is_recovering = np.full(len(fft_data), False, dtype=bool)
        for tstamp in discharge_timestamps:
            is_recovering = is_recovering | np.where(
                (
                    ((fft_data["timestamp"] - tstamp) < 0.01)
                    & ((fft_data["timestamp"] - tstamp) > 0)
                ),
                True,
                False,
            )
        fft_data["is_recovering"] = is_recovering

        hit_dict_fft = {}
        plot_dict_fft = {}
        cut_data = fft_data.query("is_recovering==0")
        log.debug(f"cut_data shape: {len(cut_data)}")
        for name, cut in kwarg_dict_fft["cut_parameters"].items():
            cut_dict, cut_plots = generate_cut_classifiers(
                cut_data,
                {name: cut},
                kwarg_dict.get("rounding", 4),
                display=1 if args.plot_path else 0,
            )
            hit_dict_fft.update(cut_dict)
            plot_dict_fft.update(cut_plots)

            log.debug(f"{name} calculated cut_dict is: {json.dumps(cut_dict, indent=2)}")

            ct_mask = np.full(len(cut_data), True, dtype=bool)
            for outname, info in cut_dict.items():
                # convert to pandas eval
                exp = info["expression"]
                for key in info.get("parameters", None):
                    exp = re.sub(f"(?<![a-zA-Z0-9]){key}(?![a-zA-Z0-9])", f"@{key}", exp)
                cut_data[outname] = cut_data.eval(exp, local_dict=info.get("parameters", None))
                if "_classifier" not in outname:
                    ct_mask = ct_mask & cut_data[outname]
            cut_data = cut_data[ct_mask]

        log.debug("fft cuts applied")
        log.debug(f"cut_dict is: {json.dumps(hit_dict_fft, indent=2)}")

    else:
        hit_dict_fft = {}
        plot_dict_fft = {}

    kwarg_dict_cal = kwarg_dict["cal_fields"]

    cut_fields = get_keys(
        [
            key.replace(f"{args.channel}/dsp/", "")
            for key in ls(args.cal_files[0], f"{args.channel}/dsp/")
        ],
        kwarg_dict_cal["cut_parameters"],
    )
    if "initial_cal_cuts" in kwarg_dict:
        init_cal = kwarg_dict["initial_cal_cuts"]
        cut_fields += get_keys(
            [
                key.replace(f"{args.channel}/dsp/", "")
                for key in ls(args.cal_files[0], f"{args.channel}/dsp/")
            ],
            init_cal["cut_parameters"],
        )

    # load data in
    data, threshold_mask = load_data(
        args.cal_files,
        f"{args.channel}/dsp",
        {},
        [*cut_fields, "timestamp", "trapTmax"],
        threshold=kwarg_dict_cal.get("threshold", 0),
        return_selection_mask=True,
        cal_energy_param="trapTmax",
    )

    if args.pulser_file:
        pulser_dict = Props.read_from(args.pulser_file)
        mask = np.array(pulser_dict["mask"])

    elif args.tcm_filelist:
        # get pulser mask from tcm files
        with open(args.tcm_filelist) as f:
            tcm_files = f.read().splitlines()
        tcm_files = sorted(np.unique(tcm_files))
        ids, mask = get_tcm_pulser_ids(
            tcm_files, args.channel, kwarg_dict["pulser_multiplicity_threshold"]
        )
    else:
        msg = "No pulser file or tcm filelist provided"
        raise ValueError(msg)

    data["is_pulser"] = mask[threshold_mask]

    discharges = data["t_sat_lo"] > 0
    discharge_timestamps = np.where(data["timestamp"][discharges])[0]
    is_recovering = np.full(len(data), False, dtype=bool)
    for tstamp in discharge_timestamps:
        is_recovering = is_recovering | np.where(
            (((data["timestamp"] - tstamp) < 0.01) & ((data["timestamp"] - tstamp) > 0)),
            True,
            False,
        )
    data["is_recovering"] = is_recovering

    rng = np.random.default_rng()
    mask = np.full(len(data.query("~is_pulser & ~is_recovering")), False, dtype=bool)
    mask[rng.choice(len(data.query("~is_pulser & ~is_recovering")), 4000, replace=False)] = True

    if "initial_cal_cuts" in kwarg_dict:
        init_cal = kwarg_dict["initial_cal_cuts"]
        hit_dict_init_cal, plot_dict_init_cal = generate_cut_classifiers(
            data.query("~is_pulser & ~is_recovering")[mask],
            init_cal["cut_parameters"],
            init_cal.get("rounding", 4),
            display=1 if args.plot_path else 0,
        )
        ct_mask = np.full(len(data), True, dtype=bool)
        for outname, info in hit_dict_init_cal.items():
            # convert to pandas eval
            exp = info["expression"]
            for key in info.get("parameters", None):
                exp = re.sub(f"(?<![a-zA-Z0-9]){key}(?![a-zA-Z0-9])", f"@{key}", exp)
            data[outname] = data.eval(exp, local_dict=info.get("parameters", None))
            ct_mask = ct_mask & data[outname]

        data = data[ct_mask]
        mask = mask[ct_mask]

    else:
        hit_dict_init_cal = {}
        plot_dict_init_cal = {}

    if len(data.query("is_pulser & ~is_recovering")) > 500:
        data = data.query("is_pulser & ~is_recovering")
    else:
        data = data.query("~is_pulser & ~is_recovering")[mask]

    hit_dict_cal, plot_dict_cal = generate_cut_classifiers(
        data,
        kwarg_dict_cal["cut_parameters"],
        kwarg_dict.get("rounding", 4),
        display=1 if args.plot_path else 0,
    )

    hit_dict = {**hit_dict_fft, **hit_dict_init_cal, **hit_dict_cal}
    plot_dict = {**plot_dict_fft, **plot_dict_init_cal, **plot_dict_cal}

    pathlib.Path(os.path.dirname(args.save_path)).mkdir(parents=True, exist_ok=True)
    Props.write_to(args.save_path, hit_dict)

    if args.plot_path:
        pathlib.Path(os.path.dirname(args.plot_path)).mkdir(parents=True, exist_ok=True)
        with open(args.plot_path, "wb") as f:
            pkl.dump({"qc": plot_dict}, f, protocol=pkl.HIGHEST_PROTOCOL)
