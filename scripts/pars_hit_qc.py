from __future__ import annotations

import argparse
import json
import logging
import pickle as pkl
import re
import warnings
from pathlib import Path

import numpy as np
from legendmeta import LegendMetadata, TextDB
from legendmeta.catalog import Props
from lgdo.lh5 import ls
from pygama.pargen.data_cleaning import (
    generate_cut_classifiers,
    get_keys,
    get_tcm_pulser_ids,
)
from pygama.pargen.utils import load_data
from util.convert_np import convert_dict_np_to_float

log = logging.getLogger(__name__)

warnings.filterwarnings(action="ignore", category=RuntimeWarning)


if __name__ == "__main__":
    argparser = argparse.ArgumentParser()
    argparser.add_argument("--cal_files", help="cal_files", nargs="*", type=str)
    argparser.add_argument("--fft_files", help="fft_files", nargs="*", type=str)

    argparser.add_argument("--tcm_filelist", help="tcm_filelist", type=str, required=False)
    argparser.add_argument("--pulser_file", help="pulser_file", type=str, required=False)
    argparser.add_argument(
        "--overwrite_files",
        help="overwrite_files",
        type=str,
        required=False,
        nargs="*",
    )

    argparser.add_argument("--configs", help="config", type=str, required=True)
    argparser.add_argument("--metadata", help="metadata", type=str, required=True)
    argparser.add_argument("--log", help="log_file", type=str)

    argparser.add_argument("--datatype", help="Datatype", type=str, required=True)
    argparser.add_argument("--timestamp", help="Timestamp", type=str, required=True)
    argparser.add_argument("--channel", help="Channel", type=str, required=True)
    argparser.add_argument("--tier", help="tier", type=str, default="hit")

    argparser.add_argument("--plot_path", help="plot_path", type=str, required=False)
    argparser.add_argument("--save_path", help="save_path", type=str)
    args = argparser.parse_args()

    configs = TextDB(args.configs, lazy=True).on(args.timestamp, system=args.datatype)
    config_dict = configs["snakemake_rules"]["pars_hit_qc"]
    if "logging" in config_dict["options"]:
        log_config = config_dict["options"]["logging"]
        log_config = Props.read_from(log_config)
        if args.log is not None:
            Path(args.log).parent.mkdir(parents=True, exist_ok=True)
            log_config["handlers"]["file"]["filename"] = args.log
        logging.config.dictConfig(log_config)
        log = logging.getLogger(config_dict["options"].get("logger", "prod"))
    else:
        if args.log is not None:
            Path(args.log).parent.makedir(parents=True, exist_ok=True)
            logging.basicConfig(level=logging.INFO, filename=args.log, filemode="w")
        log = logging.getLogger(__name__)

    meta = LegendMetadata(path=args.metadata)
    chmap = meta.channelmap(args.timestamp, system=args.datatype)
    channel = f"ch{chmap[args.channel].daq.rawid:07}"

    # get metadata dictionary
    channel_dict = config_dict["inputs"]["qc_config"][args.channel]
    kwarg_dict = Props.read_from(channel_dict)

    if args.overwrite_files:
        overwrite = Props.read_from(args.overwrite_files)
        if channel in overwrite:
            overwrite = overwrite[channel]["pars"]["operations"]
        else:
            overwrite = None
    else:
        overwrite = None

    if len(args.fft_files) == 1 and Path(args.fft_files[0]).suffix == ".filelist":
        with Path(args.fft_files[0]).open() as f:
            fft_files = f.read().splitlines()
    else:
        fft_files = args.fft_files

    if len(args.cal_files) == 1 and Path(args.cal_files[0]).suffix == ".filelist":
        with Path(args.cal_files[0]).open() as f:
            cal_files = f.read().splitlines()
    else:
        cal_files = args.fft_files

    kwarg_dict_fft = kwarg_dict["fft_fields"]
    if len(fft_files) > 0:
        fft_fields = get_keys(
            [key.replace(f"{channel}/dsp/", "") for key in ls(fft_files[0], f"{channel}/dsp/")],
            kwarg_dict_fft["cut_parameters"],
        )

        fft_data = load_data(
            fft_files,
            f"{channel}/dsp",
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

    if overwrite is not None:
        for name in kwarg_dict_fft["cut_parameters"]:
            for cut_name, cut_dict in overwrite.items():
                if name in cut_name:
                    hit_dict_fft.update({cut_name: cut_dict})

    kwarg_dict_cal = kwarg_dict["cal_fields"]

    cut_fields = get_keys(
        [key.replace(f"{channel}/dsp/", "") for key in ls(cal_files[0], f"{channel}/dsp/")],
        kwarg_dict_cal["cut_parameters"],
    )
    if "initial_cal_cuts" in kwarg_dict:
        init_cal = kwarg_dict["initial_cal_cuts"]
        cut_fields += get_keys(
            [key.replace(f"{channel}/dsp/", "") for key in ls(cal_files[0], f"{channel}/dsp/")],
            init_cal["cut_parameters"],
        )

    # load data in
    data, threshold_mask = load_data(
        cal_files,
        f"{channel}/dsp",
        {},
        [*cut_fields, "timestamp", "trapTmax", "t_sat_lo"],
        threshold=kwarg_dict_cal.get("threshold", 0),
        return_selection_mask=True,
        cal_energy_param="trapTmax",
    )

    if args.pulser_file:
        pulser_dict = Props.read_from(args.pulser_file)
        mask = np.array(pulser_dict["mask"])

    elif args.tcm_filelist:
        # get pulser mask from tcm files
        with Path(args.tcm_filelist).open() as f:
            tcm_files = f.read().splitlines()
        tcm_files = sorted(np.unique(tcm_files))
        ids, mask = get_tcm_pulser_ids(
            tcm_files, channel, kwarg_dict["pulser_multiplicity_threshold"]
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
            if "classifier" not in outname:
                ct_mask = ct_mask & data[outname]

        mask = mask[ct_mask[(~data["is_pulser"] & ~data["is_recovering"]).to_numpy()]]
        data = data[ct_mask]
        log.debug("initial cal cuts applied")
        log.debug(f"cut_dict is: {json.dumps(hit_dict_init_cal, indent=2)}")

    else:
        hit_dict_init_cal = {}
        plot_dict_init_cal = {}

    if len(data.query("is_pulser & ~is_recovering")) < 500:
        data = data.query("is_pulser & ~is_recovering")
    else:
        data = data.query("~is_pulser & ~is_recovering")[mask]

    hit_dict_cal, plot_dict_cal = generate_cut_classifiers(
        data,
        kwarg_dict_cal["cut_parameters"],
        kwarg_dict.get("rounding", 4),
        display=1 if args.plot_path else 0,
    )

    if overwrite is not None:
        for name in kwarg_dict_cal["cut_parameters"]:
            for cut_name, cut_dict in overwrite.items():
                if name in cut_name:
                    hit_dict_cal.update({cut_name: cut_dict})

    hit_dict = {**hit_dict_fft, **hit_dict_init_cal, **hit_dict_cal}
    plot_dict = {**plot_dict_fft, **plot_dict_init_cal, **plot_dict_cal}

    hit_dict = convert_dict_np_to_float(hit_dict)

    Path(args.save_path).parent.mkdir(parents=True, exist_ok=True)
    Props.write_to(args.save_path, hit_dict)

    if args.plot_path:
        Path(args.plot_path).parent.mkdir(parents=True, exist_ok=True)
        with Path(args.plot_path).open("wb") as f:
            pkl.dump({"qc": plot_dict}, f, protocol=pkl.HIGHEST_PROTOCOL)
