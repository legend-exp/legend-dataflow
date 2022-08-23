import numpy as np
import os,json
import pathlib
import argparse
import logging

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


files = sorted(args.files)

with open(args.ecal_file, 'r') as o:
    cal_dict = json.load(o)

with open(args.eres_file, 'r') as o:
    eres_dict = json.load(o)

energy_param = 'cuspEmax_cal'
cal_energy_param = 'cuspEmax_ctc_cal'

eres_pars = [eres_dict['cuspEmax_ctc_cal']['m0'],eres_dict['cuspEmax_ctc_cal']['m1']]
    

cal_dict, out_dict = cal_aoe(files, f'{args.channel}/dsp',cal_dict, energy_param, cal_energy_param, eres_pars,
                            dt_corr=False, cut_parameters={}, plot_savepath=args.plot_file)

outputs= [ "cuspEmax_ctc_cal", "zacEmax_ctc_cal", "trapEmax_ctc_cal", 
            "AoE_Classifier", "AoE_Low_Cut", "AoE_Double_Sided_Cut", "Quality_cuts"]
final_hit_dict = {"outputs":outputs, "operations":cal_dict}


pathlib.Path(os.path.dirname(args.hit_pars)).mkdir(parents=True, exist_ok=True)
with open(args.hit_pars, 'w') as w:
    json.dump(final_hit_dict,w, indent=4)

pathlib.Path(os.path.dirname(args.aoe_results)).mkdir(parents=True, exist_ok=True)
with open(args.aoe_results, 'w') as w:
    json.dump(out_dict,w, indent=4)