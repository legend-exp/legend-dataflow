from __future__ import annotations

import argparse
import json
import logging
import os
import pathlib
import pickle as pkl
import re
import warnings

import lgdo.lh5 as lh5
import numpy as np
from legendmeta import LegendMetadata
from legendmeta.catalog import Props
from lgdo.lh5 import ls
from pygama.pargen.data_cleaning import (
    generate_cut_classifiers,
    get_keys,
)

log = logging.getLogger(__name__)

warnings.filterwarnings(action="ignore", category=RuntimeWarning)


if __name__ == "__main__":
    argparser = argparse.ArgumentParser()
    argparser.add_argument("--phy_files", help="cal_files", nargs="*", type=str)

    argparser.add_argument("--configs", help="config", type=str, required=True)
    argparser.add_argument("--datatype", help="Datatype", type=str, required=True)
    argparser.add_argument("--timestamp", help="Timestamp", type=str, required=True)
    argparser.add_argument("--channel", help="Channel", type=str, required=True)

    argparser.add_argument("--log", help="log_file", type=str)

    argparser.add_argument("--plot_path", help="plot_path", type=str, nargs="*", required=False)
    argparser.add_argument(
        "--save_path",
        help="save_path",
        type=str,
        nargs="*",
    )
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
    channel_dict = channel_dict["pars_pht_qc"]["inputs"]["qc_config"][args.channel]

    sto = lh5.LH5Store()

    # sort files in dictionary where keys are first timestamp from run
    bl_mask = np.array([], dtype=bool)
    if isinstance(args.phy_files, list):
        phy_files = []
        for file in sorted(args.phy_files):
            with open(file) as f:
                run_files = f.read().splitlines()
            if len(run_files) == 0:
                continue
            else:
                run_files = sorted(np.unique(run_files))
                phy_files += run_files
                bls = sto.read("ch1027200/dsp/", run_files, field_mask=["wf_max", "bl_mean"])[0]
                puls = sto.read("ch1027201/dsp/", run_files, field_mask=["trapTmax"])[0]
                bl_idxs = ((bls["wf_max"].nda - bls["bl_mean"].nda) > 1000) & (
                    puls["trapTmax"].nda < 200
                )
                bl_mask = np.append(bl_mask, bl_idxs)
    else:
        with open(args.phy_files) as f:
            phy_files = f.read().splitlines()
        phy_files = sorted(np.unique(phy_files))
        bls = sto.read("ch1027200/dsp/", phy_files, field_mask=["wf_max", "bl_mean"])[0]
        puls = sto.read("ch1027201/dsp/", phy_files, field_mask=["trapTmax"])[0]
        bl_mask = ((bls["wf_max"].nda - bls["bl_mean"].nda) > 1000) & (puls["trapTmax"].nda < 200)

    kwarg_dict = Props.read_from(channel_dict)
    kwarg_dict_fft = kwarg_dict["fft_fields"]

    cut_fields = get_keys(
        [
            key.replace(f"{args.channel}/dsp/", "")
            for key in ls(phy_files[0], f"{args.channel}/dsp/")
        ],
        kwarg_dict_fft["cut_parameters"],
    )

    data = sto.read(
        f"{args.channel}/dsp/",
        phy_files,
        field_mask=[*cut_fields, "daqenergy", "t_sat_lo", "timestamp"],
        idx=np.where(bl_mask)[0],
    )[0].view_as("pd")

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

    log.debug(f"{len(discharge_timestamps)} discharges found in {len(data)} events")

    hit_dict = {}
    plot_dict = {}
    cut_data = data.query("is_recovering==0")
    log.debug(f"cut_data shape: {len(cut_data)}")
    for name, cut in kwarg_dict_fft["cut_parameters"].items():
        cut_dict, cut_plots = generate_cut_classifiers(
            cut_data,
            {name: cut},
            kwarg_dict.get("rounding", 4),
            display=1 if args.plot_path else 0,
        )
        hit_dict.update(cut_dict)
        plot_dict.update(cut_plots)

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
    log.debug(f"cut_dict is: {json.dumps(hit_dict, indent=2)}")

    for file in args.save_path:
        pathlib.Path(os.path.dirname(file)).mkdir(parents=True, exist_ok=True)
        Props.write_to(file, {"pars": {"operations": hit_dict}})

    if args.plot_path:
        for file in args.plot_path:
            pathlib.Path(os.path.dirname(file)).mkdir(parents=True, exist_ok=True)
            with open(file, "wb") as f:
                pkl.dump({"qc": plot_dict}, f, protocol=pkl.HIGHEST_PROTOCOL)
