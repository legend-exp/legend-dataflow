import argparse, os, pathlib
from pygama.dsp.utils import numba_defaults

numba_defaults.cache = False
numba_defaults.boundscheck = True

from pygama.dsp.build_dsp import build_dsp

import numpy as np

import json
from collections import OrderedDict
from util.metadata_loading import *
import logging

argparser = argparse.ArgumentParser()
argparser.add_argument("--configs", help="configs path", type=str, required=True)
argparser.add_argument("--datatype", help="Datatype", type=str, required=True)
argparser.add_argument("--timestamp", help="Timestamp", type=str, required=True)
argparser.add_argument("--pars_file", help="database file for detector", type=str)
argparser.add_argument("--log", help="log file", type=str)
argparser.add_argument("input", help="input file", type=str)
argparser.add_argument("output", help="output file", type=str)
args = argparser.parse_args()

logging.basicConfig(level=logging.DEBUG, filename=args.log, filemode='w')
logging.getLogger('numba').setLevel(logging.INFO)
logging.getLogger('parse').setLevel(logging.INFO)

cfg_file = os.path.join(args.configs, 'key_resolve.jsonl')

channel_dict = config_catalog.get_config(cfg_file, args.configs, args.timestamp, args.datatype)

channel_dict = channel_dict['snakemake_rules']['tier_dsp']["inputs"]['processing_chain']

with open(args.pars_file, 'r') as db:
    database_dic = json.load(db)

pathlib.Path(os.path.dirname(args.output)).mkdir(parents=True, exist_ok=True)

rand_num = f'{np.random.randint(0,99999):05d}'
temp_output = f'{args.output}.{rand_num}'



build_dsp(args.input, temp_output, {}, database = database_dic, chan_config=channel_dict, write_mode='r')

os.rename(temp_output, args.output)
