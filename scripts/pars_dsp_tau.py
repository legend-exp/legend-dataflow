from pygama.dsp.utils import numba_defaults

numba_defaults.cache = False
numba_defaults.boundscheck = True

import argparse, os, pathlib
import pygama
import pygama.pargen.extract_tau as dpp
from legendmeta import LegendMetadata
import logging
import pickle as pkl

import json
from collections import OrderedDict 

argparser = argparse.ArgumentParser()
argparser.add_argument("--configs", help="configs path", type=str, required=True)
argparser.add_argument("--log", help="log file", type=str)
argparser.add_argument("--datatype", help="Datatype", type=str, required=True)
argparser.add_argument("--timestamp", help="Timestamp", type=str, required=True)
argparser.add_argument("--channel", help="Channel", type=str, required=True)
argparser.add_argument("--plot_path", help="plot path", type=str, required=False)
argparser.add_argument("--output_file", help="output file", type=str, required=True)
argparser.add_argument("input", help="input files", nargs='*',type=str)
args = argparser.parse_args()

logging.basicConfig(level=logging.DEBUG, filename=args.log, filemode='w')
logging.getLogger('numba').setLevel(logging.INFO)
logging.getLogger('parse').setLevel(logging.INFO)

configs = LegendMetadata(path = args.configs)
config_dict = configs.on(args.timestamp, system=args.datatype)
channel_dict = config_dict['snakemake_rules']['pars_dsp_tau']["inputs"]['processing_chain'][args.channel]
kwarg_dict = config_dict['snakemake_rules']['pars_dsp_tau']["inputs"]['tau_config'][args.channel] 

with open(kwarg_dict,"r") as r:
    kwarg_dict = json.load(r)

if kwarg_dict["run_tau"]==True:
    kwarg_dict.pop("run_tau")
    input_file = args.input
    if isinstance(input_file, list):
        if input_file[0].split('.')[-1] == 'filelist':
            input_file = args.input[0]
            with open(input_file) as f:
                input_file = f.read().splitlines()

    
    if args.plot_path:
        out_dict,plot_dict = dpp.dsp_preprocess_decay_const(input_file, channel_dict, f'{args.channel}/raw', **kwarg_dict, display=1) 
        pathlib.Path(os.path.dirname(args.plot_path)).mkdir(parents=True, exist_ok=True)
        with open(args.plot_path,"wb") as f:
            pkl.dump(plot_dict,f)
    else:
        out_dict = dpp.dsp_preprocess_decay_const(input_file, channel_dict, f'{args.channel}/raw', **kwarg_dict) 
else:
    out_dict = {}

pathlib.Path(os.path.dirname(args.output_file)).mkdir(parents=True, exist_ok=True)
with open(args.output_file,"w") as f:
    json.dump(out_dict,f, indent=4)