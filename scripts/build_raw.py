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
argparser.add_argument("--log", help="log file", type=str)
args = argparser.parse_args()

logging.basicConfig(level=logging.INFO, filename=args.log, filemode='w')

pathlib.Path(os.path.dirname(args.output)).mkdir(parents=True, exist_ok=True)

configs = LegendMetadata(path = args.configs)
channel_dict = configs.on(args.timestamp, system=args.datatype)['snakemake_rules']['tier_raw']["inputs"]['raw_config']
with open(channel_dict, "r") as f:
  out_spec = json.load(f)

rand_num = f'{np.random.randint(0,99999):05d}'
temp_output = f'{args.output}.{rand_num}'

build_raw(args.input, in_stream_type='ORCA', out_spec=out_spec, filekey = temp_output, buffer_size=1024)

os.rename(temp_output, args.output)
