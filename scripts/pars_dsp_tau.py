import argparse, os, pathlib
import pygama
import pygama.pargen.extract_tau as dpp
from util.metadata_loading import *

import json
from collections import OrderedDict 

argparser = argparse.ArgumentParser()
argparser.add_argument("--configs", help="configs path", type=str, required=True)
argparser.add_argument("--datatype", help="Datatype", type=str, required=True)
argparser.add_argument("--timestamp", help="Timestamp", type=str, required=True)
argparser.add_argument("--channel", help="Channel", type=str, required=True)
argparser.add_argument("--plot_path", help="plot path", type=str, required=False)
argparser.add_argument("--output_file", help="output file", type=str, required=True)

argparser.add_argument("input", help="input files", nargs='*',type=str)
args = argparser.parse_args()

cfg_file = os.path.join(args.configs, 'key_resolve.jsonl')
channel_dict = config_catalog.get_config(cfg_file, args.configs, args.timestamp, args.datatype)
channel_dict = channel_dict['snakemake_rules']['pars_dsp_tau']["inputs"]['processing_chain'][args.channel]

print(channel_dict)

with open(channel_dict , 'r') as config_file:
    config_dic = json.load(config_file, object_pairs_hook=OrderedDict)

input_file = args.input
if isinstance(input_file, list):
    if input_file[0].split('.')[-1] == 'filelist':
        input_file = args.input[0]
        with open(input_file) as f:
            input_file = f.read().splitlines()


print(input_file)

pathlib.Path(os.path.dirname(args.output_file)).mkdir(parents=True, exist_ok=True)
if args.plot_path:
    out_dict = dpp.dsp_preprocess_decay_const(input_file, config_dic, args.channel, plot_path=args.plot_path, verbose=True) 
else:
    out_dict = dpp.dsp_preprocess_decay_const(input_file, config_dic, args.channel,verbose=True) 

with open(args.output_file,"w") as f:
    json.dump(out_dict,f, indent=4)