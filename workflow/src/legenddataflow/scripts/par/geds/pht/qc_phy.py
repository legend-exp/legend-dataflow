from __future__ import annotations

import argparse
import json
import pickle as pkl
import re
import warnings
from pathlib import Path

import lgdo.lh5 as lh5
import numpy as np
from dbetto import TextDB
from dbetto.catalog import Props
from lgdo.lh5 import ls
from pygama.pargen.data_cleaning import (
    generate_cut_classifiers,
    get_keys,
)

from .....convert_np import convert_dict_np_to_float
from .....log import build_log

warnings.filterwarnings(action="ignore", category=RuntimeWarning)


def par_geds_pht_qc_phy() -> None:
    argparser = argparse.ArgumentParser()
    argparser.add_argument("--phy-files", help="cal_files", nargs="*", type=str)

    argparser.add_argument("--configs", help="config", type=str, required=True)
    argparser.add_argument("--metadata", help="metadata path", type=str, required=True)
    argparser.add_argument("--log", help="log_file", type=str)

    argparser.add_argument("--datatype", help="Datatype", type=str, required=True)
    argparser.add_argument("--timestamp", help="Timestamp", type=str, required=True)
    argparser.add_argument("--channel", help="Channel", type=str, required=True)
    argparser.add_argument("--table-name", help="table name", type=str, required=True)
    argparser.add_argument(
        "--baseline-name", help="table name", type=str, required=True
    )
    argparser.add_argument("--pulser-name", help="table name", type=str, required=True)

    argparser.add_argument(
        "--plot-path", help="plot_path", type=str, nargs="*", required=False
    )
    argparser.add_argument(
        "--save-path",
        help="save_path",
        type=str,
        nargs="*",
    )
    args = argparser.parse_args()

    configs = TextDB(args.configs, lazy=True).on(args.timestamp, system=args.datatype)
    config_dict = configs["snakemake_rules"]["pars_pht_qc"]

    log = build_log(config_dict, args.log)

    # get metadata dictionary
    channel_dict = config_dict["qc_config"][args.channel]
    kwarg_dict = Props.read_from(channel_dict)

    sto = lh5.LH5Store()

    # sort files in dictionary where keys are first timestamp from run
    bl_mask = np.array([], dtype=bool)
    if isinstance(args.phy_files, list):
        phy_files = []
        for file in sorted(args.phy_files):
            with Path(file).open() as f:
                run_files = f.read().splitlines()
            if len(run_files) == 0:
                continue
            run_files = sorted(np.unique(run_files))
            phy_files += run_files
            bls = sto.read(
                args.baseline_name, run_files, field_mask=["wf_max", "bl_mean"]
            )[0]
            puls = sto.read(args.pulser_name, run_files, field_mask=["trapTmax"])[0]
            bl_idxs = ((bls["wf_max"].nda - bls["bl_mean"].nda) > 1000) & (
                puls["trapTmax"].nda < 200
            )
            bl_mask = np.append(bl_mask, bl_idxs)
    else:
        with Path(args.phy_files).open() as f:
            phy_files = f.read().splitlines()
        phy_files = sorted(np.unique(phy_files))
        bls = sto.read(args.baseline_name, phy_files, field_mask=["wf_max", "bl_mean"])[
            0
        ]
        puls = sto.read(args.pulser_name, phy_files, field_mask=["trapTmax"])[0]
        bl_mask = ((bls["wf_max"].nda - bls["bl_mean"].nda) > 1000) & (
            puls["trapTmax"].nda < 200
        )

    kwarg_dict_fft = kwarg_dict["fft_fields"]

    search_name = (
        args.table_name if args.table_name[-1] == "/" else args.table_name + "/"
    )

    cut_fields = get_keys(
        [key.replace(search_name, "") for key in ls(phy_files[0], search_name)],
        kwarg_dict_fft["cut_parameters"],
    )

    data = sto.read(
        args.table_name,
        phy_files,
        field_mask=[*cut_fields, "daqenergy", "t_sat_lo", "timestamp"],
        idx=np.where(bl_mask)[0],
    )[0].view_as("pd")

    discharges = data["t_sat_lo"] > 0
    discharge_timestamps = np.where(data["timestamp"][discharges])[0]
    is_recovering = np.full(len(data), False, dtype=bool)
    for tstamp in discharge_timestamps:
        is_recovering = is_recovering | np.where(
            (
                ((data["timestamp"] - tstamp) < 0.01)
                & ((data["timestamp"] - tstamp) > 0)
            ),
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
            cut_data[outname] = cut_data.eval(
                exp, local_dict=info.get("parameters", None)
            )
            if "_classifier" not in outname:
                ct_mask = ct_mask & cut_data[outname]
        cut_data = cut_data[ct_mask]

    log.debug("fft cuts applied")
    log.debug(f"cut_dict is: {json.dumps(hit_dict, indent=2)}")

    hit_dict = convert_dict_np_to_float(hit_dict)

    for file in args.save_path:
        Path(file).name.mkdir(parents=True, exist_ok=True)
        Props.write_to(file, {"pars": {"operations": hit_dict}})

    if args.plot_path:
        for file in args.plot_path:
            Path(file).parent.mkdir(parents=True, exist_ok=True)
            with Path(file).open("wb") as f:
                pkl.dump({"qc": plot_dict}, f, protocol=pkl.HIGHEST_PROTOCOL)
