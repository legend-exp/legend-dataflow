from __future__ import annotations

import argparse
import json
import logging
import os
import pathlib
import pickle as pkl
import warnings

from legendmeta import LegendMetadata
from legendmeta.catalog import Props
from pygama.pargen.utils import get_tcm_pulser_ids, load_data
from pygama.pargen.cuts import generate_cuts

log = logging.getLogger(__name__)

warnings.filterwarnings(action="ignore", category=RuntimeWarning)


if __name__ == "__main__":
    argparser = argparse.ArgumentParser()
    argparser.add_argument("--files", help="files", nargs="*", type=str)
    argparser.add_argument("--tcm_filelist", help="tcm_filelist", type=str, required=True)

    argparser.add_argument("--configs", help="config", type=str, required=True)
    argparser.add_argument("--datatype", help="Datatype", type=str, required=True)
    argparser.add_argument("--timestamp", help="Timestamp", type=str, required=True)
    argparser.add_argument("--channel", help="Channel", type=str, required=True)
    argparser.add_argument("--tier", help="tier", type=str, default="hit")

    argparser.add_argument("--log", help="log_file", type=str)

    argparser.add_argument("--plot_path", help="plot_path", type=str, required=False, nargs="*")
    argparser.add_argument("--save_path", help="save_path", type=str, nargs="*")
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
    if args.tier == "hit":
        channel_dict = channel_dict["pars_hit_qc"]["inputs"]["ecal_config"][args.channel]
    elif args.tier == "pht":
        channel_dict = channel_dict["pars_pht_qc"]["inputs"]["ecal_config"][args.channel]
    else:
        msg = "invalid tier"
        raise ValueError(msg)

    kwarg_dict = Props.read_from(channel_dict)

    # load data in
    data, threshold_mask = load_data(
        args.files,
        f"{args.channel}/dsp",
        hit_dict,
        list(kwarg_dict["cut_parameters"])
        + ["timestamp", "trapTmax"],
        threshold=kwarg_dict["threshold"],
        return_selection_mask=True,
        cal_energy_param="trapTmax",
    )

    # get pulser mask from tcm files
    with open(args.tcm_filelist) as f:
        tcm_files = f.read().splitlines()
    tcm_files = sorted(np.unique(tcm_files))
    ids, mask = get_tcm_pulser_ids(
        tcm_files, args.channel, kwarg_dict.pop("pulser_multiplicity_threshold")
    )
    data["is_pulser"] = mask[threshold_mask]

    hit_dict, plot_dict = generate_cuts(
            data,
            cut_dict,
            kwarg_dict.get("rounding",4),
            display=1 if args.plot_path else 0,
        )
    if isinstance(args.save_path, string):
        save_path = [args.save_path]
    else:
        save_path = args.save_path
    for file in save_path
        pathlib.Path(os.path.dirname(save_path)).mkdir(parents=True, exist_ok=True)
        with open(file, "w") as f:
            json.dump(hit_dict, f, indent=4)

    if args.plot_path:
        if isinstance(args.plot_path, string):
            plot_path = [args.plot_path]
        else:
            plot_path = args.plot_path
        for file in plot_path:
             pathlib.Path(os.path.dirname(plot_path)).mkdir(parents=True, exist_ok=True)
            with open(plot_path, "wb") as f:
                pkl.dump({"qc":plot_dict}, f, protocol=pkl.HIGHEST_PROTOCOL)