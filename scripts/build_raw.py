import argparse, os, pathlib, json
import logging

from util.metadata_loading import *
from util.Props import *

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

# ToDo: Atomic file creation

cfg_file = os.path.join(args.configs, 'key_resolve.jsonl')
channel_dict = config_catalog.get_config(cfg_file, args.configs, args.timestamp, args.datatype)
channel_dict = channel_dict['snakemake_rules']['tier_raw']["inputs"]['raw_config']
with open(channel_dict, "r") as f:
  out_spec = json.load(f)

rand_num = f'{np.random.randint(0,99999):05d}'
temp_output = f'{args.output}.{rand_num}'

build_raw(args.input, in_stream_type='ORCA', out_spec=json.loads(out_spec), filekey = temp_output, buffer_size=1024)

os.rename(temp_output, args.output)
