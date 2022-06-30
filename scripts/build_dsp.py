import argparse, os, pathlib

import pygama
from pygama.dsp.build_dsp import build_dsp

import json
from collections import OrderedDict
from util.metadata_loading import *

argparser = argparse.ArgumentParser()
argparser.add_argument("--configs", help="configs path", type=str, required=True)
argparser.add_argument("--datatype", help="Datatype", type=str, required=True)
argparser.add_argument("--timestamp", help="Timestamp", type=str, required=True)
argparser.add_argument("--pars_file", help="database file for detector", type=str)
argparser.add_argument("input", help="input file", type=str)
argparser.add_argument("output", help="output file", type=str)
args = argparser.parse_args()

cfg_file = os.path.join(args.configs, 'key_resolve.jsonl')

channel_dict = config_catalog.get_config(cfg_file, args.configs, args.timestamp, args.datatype)

channel_dict = channel_dict['snakemake_rules']['tier_dsp']["inputs"]['processing_chain']

pathlib.Path(os.path.dirname(args.output)).mkdir(parents=True, exist_ok=True)

build_dsp(args.input, args.output, {}, chan_config=channel_dict, verbose=True, overwrite=False) #database = database_dic,