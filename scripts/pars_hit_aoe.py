import numpy as np
import os,json
import pathlib
import argparse
import logging

from util.metadata_loading import *

from pygama.pargen.AoE_cal import cal_aoe


argparser = argparse.ArgumentParser()
argparser.add_argument("files", help="files", nargs='*',type=str)
argparser.add_argument("--ecal_file", help="ecal_file",type=str, required=True)
argparser.add_argument("--eres_file", help="eres_file",type=str, required=True)

argparser.add_argument("--configs", help="configs", type=str, required=True)
argparser.add_argument("--datatype", help="Datatype", type=str, required=True)
argparser.add_argument("--timestamp", help="Timestamp", type=str, required=True)
argparser.add_argument("--channel", help="Channel", type=str, required=True)

argparser.add_argument("--log", help="log_file", type=str)

argparser.add_argument("--plot_file", help="plot_file",type=str)
argparser.add_argument("--hit_pars", help="hit_pars",type=str)
argparser.add_argument("--aoe_results", help="aoe_results",type=str)
args = argparser.parse_args()

logging.basicConfig(level=logging.DEBUG, filename=args.log, filemode='w')
logging.getLogger('numba').setLevel(logging.INFO)
logging.getLogger('parse').setLevel(logging.INFO)

with open(args.files[0]) as f:
    files = f.read().splitlines()
files = sorted(files)

with open(args.ecal_file, 'r') as o:
    cal_dict = json.load(o)

cfg_file = os.path.join(args.configs, 'key_resolve.jsonl')
channel_dict = config_catalog.get_config(cfg_file, args.configs, args.timestamp, args.datatype)
channel_dict = channel_dict['snakemake_rules']['pars_hit_aoecal']["inputs"]['aoecal_config'][args.channel]

with open(channel_dict,"r") as r:
    kwarg_dict = json.load(r)

with open(args.eres_file, 'r') as o:
    eres_dict = json.load(o)

eres_pars = eres_dict[kwarg_dict["cal_energy_param"]]['eres_pars']
    

cal_dict, out_dict = cal_aoe(files, lh5_path=f'{args.channel}/dsp', cal_dict=cal_dict, eres_pars=eres_pars,
                             plot_savepath=args.plot_file, **kwarg_dict) 

#outputs= [ "cuspEmax_ctc_cal", "zacEmax_ctc_cal", "trapEmax_ctc_cal", "AoE_Corrected",
#            "AoE_Classifier", "AoE_Low_Cut", "AoE_Double_Sided_Cut", "Quality_cuts"]
final_hit_dict = { "operations":cal_dict} #"outputs":outputs,


pathlib.Path(os.path.dirname(args.hit_pars)).mkdir(parents=True, exist_ok=True)
with open(args.hit_pars, 'w') as w:
    json.dump(final_hit_dict,w, indent=4)

pathlib.Path(os.path.dirname(args.aoe_results)).mkdir(parents=True, exist_ok=True)
with open(args.aoe_results, 'w') as w:
    json.dump(out_dict,w, indent=4)