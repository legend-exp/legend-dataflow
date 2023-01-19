import argparse, os, pathlib, json
import logging

from legendmeta import LegendMetadata
from legendmeta.catalog import Props

import pygama
from pygama.raw.build_raw import * 
import numpy as np


argparser = argparse.ArgumentParser()
argparser.add_argument("input", help="input file", type=str)
argparser.add_argument("output", help="output file", type=str)
argparser.add_argument("--datatype", help="Datatype", type=str, required=True)
argparser.add_argument("--timestamp", help="Timestamp", type=str, required=True)
argparser.add_argument("--configs", help="config file", type=str)
argparser.add_argument("--chan_maps", help="chan map", type=str)
argparser.add_argument("--log", help="log file", type=str)
args = argparser.parse_args()

logging.basicConfig(level=logging.INFO, filename=args.log, filemode='w')

pathlib.Path(os.path.dirname(args.output)).mkdir(parents=True, exist_ok=True)

configs = LegendMetadata(path = args.configs)
channel_dict = configs.on(args.timestamp, system=args.datatype)['snakemake_rules']['tier_raw']["inputs"]
with open(channel_dict["gen_config"], "r") as f:
  all_config = json.load(f)  


if "ged_config" in list(channel_dict) or "spm_config" in list(channel_dict):
  with open(channel_dict["ged_config"], "r") as f:
    ged_config = json.load(f)
  with open(channel_dict["spm_config"], "r") as f:
    spm_config = json.load(f)
  hardware_configs = LegendMetadata(path = args.chan_maps)
  ge_chans = []
  spm_chans = []
  chan_map = hardware_configs.channelmaps.on(args.timestamp)
  for field in chan_map:
      if chan_map[field]["system"] == "spms":
          spm_chans.append(chan_map[field]["daq"]["fcid"])
      else:
          ge_chans.append(chan_map[field]["daq"]["fcid"])

  ged_config[list(ged_config)[0]]["geds"]["key_list"]= sorted(ge_chans)
  spm_config[list(spm_config)[0]]["spms"]["key_list"]= sorted(spm_chans)
  Props.add_to(all_config, ged_config)
  Props.add_to(all_config, spm_config)


rand_num = f'{np.random.randint(0,99999):05d}'
temp_output = f'{args.output}.{rand_num}'

build_raw(args.input, in_stream_type='ORCA', out_spec=all_config, 
          filekey = temp_output, buffer_size=1024)

os.rename(temp_output, args.output)
