import argparse
import logging
import os
import pathlib

import lgdo
import numpy as np
from daq2lh5.orca import orca_flashcam
from pygama.evt.build_tcm import *

argparser = argparse.ArgumentParser()
argparser.add_argument("input", help="input file", type=str)
argparser.add_argument("output", help="output file", type=str)
argparser.add_argument("--configs", help="config file", type=str)
argparser.add_argument("--log", help="log file", type=str)
args = argparser.parse_args()

logging.basicConfig(level=logging.DEBUG, filename=args.log, filemode="w")

pathlib.Path(os.path.dirname(args.output)).mkdir(parents=True, exist_ok=True)

rng = np.random.default_rng()
rand_num = f"{rng.integers(0,99999):05d}"
temp_output = f"{args.output}.{rand_num}"

# get the list of channels by fcid
ch_list = lgdo.ls(args.input, "/ch*")
fcid_channels = {}
for ch in ch_list:
    key = int(ch[2:])
    fcid = orca_flashcam.get_fcid(key)
    if fcid not in fcid_channels:
        fcid_channels[fcid] = []
    fcid_channels[fcid].append(f"/{ch}/raw")

# make a hardware_tcm_[fcid] for each fcid
for fcid in fcid_channels:
    out_name = f"hardware_tcm_{fcid}"
    ch_list = fcid_channels[fcid]
    build_tcm(
        [(args.input, ch_list)],
        "timestamp",
        out_file=temp_output,
        out_name=out_name,
        wo_mode="o",
    )

os.rename(temp_output, args.output)
