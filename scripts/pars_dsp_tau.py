import argparse, os, pathlib

import pygama
import pygama.pargen.get_decay_const as dpp

import json
from collections import OrderedDict 

argparser = argparse.ArgumentParser()
argparser.add_argument("--configs", help="configs path", type=str, required=True)
argparser.add_argument("--plot_path", help="plot path", type=str, required=False)
argparser.add_argument("--output_file", help="output file", type=str, required=True)

argparser.add_argument("input", help="input file", nargs='*',type=str)
args = argparser.parse_args()


f_config = os.path.join(f"{args.configs}", "main_initial_dsp_config.json")

with open(f_config) as f:
    config_dic = json.load(f, object_pairs_hook=OrderedDict) 

input_file = args.input

if isinstance(input_file, str):
    if input_file.split('.')[-1] == 'filelist':
        with open(input_file) as f:
            input_file = f.read().splitlines()[0]
elif isinstance(input_file,list):
    input_file = input_file[0]

pathlib.Path(os.path.dirname(args.output_file)).mkdir(parents=True, exist_ok=True)
#if args.plot_path:
#dpp.dsp_preprocess_decay_const(input_file, config_dic, database_file=args.output_file, plot_path=args.plot_path, verbose=True, overwrite=False) 
#else:
#dpp.dsp_preprocess_decay_const(input_file, config_dic, database_file=args.output_file, verbose=True, overwrite=False) 
out_dict = {"pz":{"tau":4}}
with open(args.output_file,"w") as f:
    json.dump(out_dict,f, indent=4)

pathlib.Path(args.plot_path).touch()