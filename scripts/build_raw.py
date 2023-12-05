import argparse
import logging
import os
import pathlib

from lgdo.utils import numba_defaults

numba_defaults.cache = False
numba_defaults.boundscheck = False

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
settings = Props.read_from(channel_dict["settings"])
channel_dict = channel_dict["out_spec"]
all_config = Props.read_from(channel_dict["gen_config"])

chmap = LegendMetadata(path=args.chan_maps)

if "geds_config" in list(channel_dict):
    ged_config = Props.read_from(channel_dict["geds_config"])

    ged_channels = list(
        chmap.channelmaps.on(args.timestamp).map("system", unique=False)["geds"].map("daq.rawid")
    )

    ged_config[list(ged_config)[0]]["geds"]["key_list"] = sorted(ged_channels)
    Props.add_to(all_config, ged_config)

if "spms_config" in list(channel_dict):
    spm_config = Props.read_from(channel_dict["spms_config"])

    spm_channels = list(
        chmap.channelmaps.on(args.timestamp).map("system", unique=False)["spms"].map("daq.rawid")
    )

    spm_config[list(spm_config)[0]]["spms"]["key_list"] = sorted(spm_channels)
    Props.add_to(all_config, spm_config)

if "auxs_config" in list(channel_dict):
    aux_config = Props.read_from(channel_dict["auxs_config"])
    aux_channels = list(
        chmap.channelmaps.on(args.timestamp).map("system", unique=False)["auxs"].map("daq.rawid")
    )
    aux_channels += list(
        chmap.channelmaps.on(args.timestamp).map("system", unique=False)["puls"].map("daq.rawid")
    )
    aux_channels += list(
        chmap.channelmaps.on(args.timestamp).map("system", unique=False)["bsln"].map("daq.rawid")
    )
    top_key = list(aux_config)[0]
    aux_config[top_key][list(aux_config[top_key])[0]]["key_list"] = sorted(aux_channels)
    Props.add_to(all_config, aux_config)

if "muon_config" in list(channel_dict):
    muon_config = Props.read_from(channel_dict["muon_config"])
    muon_channels = list(
        chmap.channelmaps.on(args.timestamp).map("system", unique=False)["muon"].map("daq.rawid")
    )
    top_key = list(muon_config)[0]
    muon_config[top_key][list(muon_config[top_key])[0]]["key_list"] = sorted(muon_channels)
    Props.add_to(all_config, muon_config)

rng = np.random.default_rng()
rand_num = f"{rng.integers(0,99999):05d}"
temp_output = f"{args.output}.{rand_num}"

build_raw(args.input, out_spec=all_config, filekey=temp_output, **settings)

os.rename(temp_output, args.output)
