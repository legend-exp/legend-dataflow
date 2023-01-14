import numpy as np
import os,json
import pathlib
import argparse
import logging
import pickle as pkl

from legendmeta import LegendMetadata

from pygama.pargen.AoE_cal import cal_aoe


argparser = argparse.ArgumentParser()
argparser.add_argument("files", help="files", nargs='*',type=str)
argparser.add_argument("--ecal_file", help="ecal_file",type=str, required=True)
argparser.add_argument("--eres_file", help="eres_file",type=str, required=True)
argparser.add_argument("--inplots", help="in_plot_path", type=str, required=False)

argparser.add_argument("--configs", help="configs", type=str, required=True)
argparser.add_argument("--datatype", help="Datatype", type=str, required=True)
argparser.add_argument("--timestamp", help="Timestamp", type=str, required=True)
argparser.add_argument("--channel", help="Channel", type=str, required=True)

argparser.add_argument("--log", help="log_file", type=str)

argparser.add_argument("--plot_file", help="plot_file",type=str, required=False)
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

configs = LegendMetadata(path = args.configs)
channel_dict = configs.on(args.timestamp, system=args.datatype)['snakemake_rules']['pars_hit_aoecal']["inputs"]['aoecal_config'][args.channel]

with open(channel_dict,"r") as r:
    kwarg_dict = json.load(r)

with open(args.eres_file, 'r') as o:
    eres_dict = json.load(o)

eres_pars = eres_dict[kwarg_dict["cal_energy_param"]]['eres_pars']

if kwarg_dict["run_aoe"] ==True:
    kwarg_dict.pop("run_aoe")
    
    if args.plot_file:
        cal_dict, out_dict,plot_dict = cal_aoe(files, lh5_path=f'{args.channel}/dsp', cal_dict=cal_dict, 
                                eres_pars=eres_pars, display=1, **kwarg_dict) 
    else:
        cal_dict, out_dict,plot_dict = cal_aoe(files, lh5_path=f'{args.channel}/dsp', cal_dict=cal_dict, 
                                eres_pars=eres_pars, **kwarg_dict) 
else:
    out_dict = {}
    plot_dict = {}

if args.plot_file:
    if args.inplots:
        with open(args.inplots, "rb") as r:
            out_plot_dict = pkl.load(r)
        out_plot_dict.update(plot_dict)
    else:
        out_plot_dict = plot_dict

    pathlib.Path(os.path.dirname(args.plot_file)).mkdir(parents=True, exist_ok=True)
    with open(args.plot_file, "wb") as w:
        pkl.dump(out_plot_dict, w)

final_hit_dict = { "operations":cal_dict}



pathlib.Path(os.path.dirname(args.hit_pars)).mkdir(parents=True, exist_ok=True)
with open(args.hit_pars, 'w') as w:
    json.dump(final_hit_dict,w, indent=4)

pathlib.Path(os.path.dirname(args.aoe_results)).mkdir(parents=True, exist_ok=True)
with open(args.aoe_results, 'w') as w:
    json.dump({"ecal": eres_dict,"aoe":out_dict},w, indent=4)