
import argparse, os, pathlib

import pygama
from pygama.io.raw_to_dsp import raw_to_dsp

import json
from collections import OrderedDict

argparser = argparse.ArgumentParser()
argparser.add_argument("--configs", help="configs path", type=str, required=True)
#argparser.add_argument("--measurement", help="Measurement", type=str, required=True)
argparser.add_argument("input", help="input file", type=str)
argparser.add_argument("--pars_file", help="database file for detector", type=str)
argparser.add_argument("output", help="output file", type=str)
args = argparser.parse_args()



f_config = os.path.join(f"{args.configs}", "main_dsp_config.json")

with open(f_config) as f:
    config_dic = json.load(f, object_pairs_hook=OrderedDict)

with open(args.pars_file) as f:
    database_dic = json.load(f, object_pairs_hook=OrderedDict)
database_dic = database_dic["1"]

pathlib.Path(os.path.dirname(args.output)).mkdir(parents=True, exist_ok=True)

#raw_to_dsp(args.input, args.output, config_dic, database = database_dic, verbose=True, overwrite=False)
pathlib.Path(args.output).touch()