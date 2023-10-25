import argparse
import logging
import os
import pathlib

import numpy as np
from daq2lh5.build_raw import build_raw
from legendmeta import LegendMetadata
from legendmeta.catalog import Props

argparser = argparse.ArgumentParser()
argparser.add_argument("input", help="input file", type=str)
argparser.add_argument("output", help="output file", type=str)
argparser.add_argument("--datatype", help="Datatype", type=str, required=True)
argparser.add_argument("--timestamp", help="Timestamp", type=str, required=True)
argparser.add_argument("--configs", help="config file", type=str)
argparser.add_argument("--chan_maps", help="chan map", type=str)
argparser.add_argument("--log", help="log file", type=str)
args = argparser.parse_args()

os.makedirs(os.path.dirname(args.log), exist_ok=True)
logging.basicConfig(level=logging.INFO, filename=args.log, filemode="w")

pathlib.Path(os.path.dirname(args.output)).mkdir(parents=True, exist_ok=True)

configs = LegendMetadata(path=args.configs)
channel_dict = configs.on(args.timestamp, system=args.datatype)["snakemake_rules"]["tier_raw"][
    "inputs"
]
all_config = Props.read_from(channel_dict["gen_config"])


if "ged_config" in list(channel_dict) or "spm_config" in list(channel_dict):
    ged_config = Props.read_from(channel_dict["ged_config"])
    spm_config = Props.read_from(channel_dict["spm_config"])

    chmap = LegendMetadata(path=args.chan_maps)
    spm_channels = list(
        chmap.channelmaps.on(args.timestamp).map("system", unique=False)["spms"].map("daq.fcid")
    )
    ged_channels = list(chmap.channelmaps.on(args.timestamp).map("daq.fcid"))
    for spm_channel in spm_channels:
        ged_channels.remove(spm_channel)

    ged_config[list(ged_config)[0]]["geds"]["key_list"] = sorted(ged_channels)
    spm_config[list(spm_config)[0]]["spms"]["key_list"] = sorted(spm_channels)
    Props.add_to(all_config, ged_config)
    Props.add_to(all_config, spm_config)

rng = np.random.default_rng()
rand_num = f"{rng.integers(0,99999):05d}"
temp_output = f"{args.output}.{rand_num}"

build_raw(
    args.input,
    in_stream_type="ORCA",
    out_spec=all_config,
    filekey=temp_output,
    buffer_size=1024,
)

os.rename(temp_output, args.output)
