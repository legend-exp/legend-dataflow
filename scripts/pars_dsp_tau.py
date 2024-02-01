import argparse
import json
import logging
import os
import pathlib
import pickle as pkl

os.environ["LGDO_CACHE"] = "false"
os.environ["LGDO_BOUNDSCHECK"] = "false"
os.environ["DSPEED_CACHE"] = "false"
os.environ["DSPEED_BOUNDSCHECK"] = "false"

import lgdo.lh5 as lh5
import numpy as np
from legendmeta import LegendMetadata
from legendmeta.catalog import Props
from pygama.pargen.extract_tau import dsp_preprocess_decay_const
from pygama.pargen.utils import get_tcm_pulser_ids

argparser = argparse.ArgumentParser()
argparser.add_argument("--configs", help="configs path", type=str, required=True)
argparser.add_argument("--log", help="log file", type=str)
argparser.add_argument("--datatype", help="Datatype", type=str, required=True)
argparser.add_argument("--timestamp", help="Timestamp", type=str, required=True)
argparser.add_argument("--channel", help="Channel", type=str, required=True)
argparser.add_argument("--plot_path", help="plot path", type=str, required=False)
argparser.add_argument("--output_file", help="output file", type=str, required=True)
argparser.add_argument("--raw_files", help="input files", nargs="*", type=str)
argparser.add_argument("--tcm_files", help="tcm_files", nargs="*", type=str)
args = argparser.parse_args()

logging.basicConfig(level=logging.DEBUG, filename=args.log, filemode="w")
logging.getLogger("numba").setLevel(logging.INFO)
logging.getLogger("parse").setLevel(logging.INFO)
logging.getLogger("lgdo").setLevel(logging.INFO)
logging.getLogger("h5py").setLevel(logging.INFO)
logging.getLogger("matplotlib").setLevel(logging.INFO)

sto = lh5.LH5Store()

configs = LegendMetadata(path=args.configs)
config_dict = configs.on(args.timestamp, system=args.datatype)
channel_dict = config_dict["snakemake_rules"]["pars_dsp_tau"]["inputs"]["processing_chain"][
    args.channel
]
kwarg_dict = config_dict["snakemake_rules"]["pars_dsp_tau"]["inputs"]["tau_config"][args.channel]

kwarg_dict = Props.read_from(kwarg_dict)

if kwarg_dict["run_tau"] is True:
    kwarg_dict.pop("run_tau")
    if isinstance(args.raw_files, list) and args.raw_files[0].split(".")[-1] == "filelist":
        input_file = args.raw_files[0]
        with open(input_file) as f:
            input_file = f.read().splitlines()
    else:
        input_file = args.raw_files

    if isinstance(args.tcm_files, list) and args.tcm_files[0].split(".")[-1] == "filelist":
        tcm_files = args.tcm_files[0]
        with open(tcm_files) as f:
            tcm_files = f.read().splitlines()
    else:
        tcm_files = args.tcm_files
    # get pulser mask from tcm files
    tcm_files = sorted(np.unique(tcm_files))
    ids, mask = get_tcm_pulser_ids(
        tcm_files, args.channel, kwarg_dict.pop("pulser_multiplicity_threshold")
    )
    data = sto.read(f"{args.channel}/raw", input_file, field_mask=["daqenergy", "timestamp"])[
        0
    ].view_as("pd")
    threshold = kwarg_dict.pop("threshold")
    cuts = np.where((data.daqenergy.to_numpy() > threshold) & (~mask))[0]

    tb_data = sto.read(
        f"{args.channel}/raw",
        input_file,
        idx=cuts,
        n_rows=kwarg_dict.pop("n_events"),
    )[0]

    out_dict, plot_dict = dsp_preprocess_decay_const(
        tb_data, channel_dict, **kwarg_dict, display=1
    )

    if args.plot_path:
        pathlib.Path(os.path.dirname(args.plot_path)).mkdir(parents=True, exist_ok=True)
        with open(args.plot_path, "wb") as f:
            pkl.dump({"tau": plot_dict}, f, protocol=pkl.HIGHEST_PROTOCOL)
else:
    out_dict = {}

pathlib.Path(os.path.dirname(args.output_file)).mkdir(parents=True, exist_ok=True)
with open(args.output_file, "w") as f:
    json.dump(out_dict, f, indent=4)
