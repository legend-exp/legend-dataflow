import numpy as np
import os,json
import pathlib
import argparse

from pygama.pargen.aoe_cal import cal_aoe


argparser = argparse.ArgumentParser()
argparser.add_argument("files", help="files", nargs='*',type=str)
argparser.add_argument("--ecal_file", help="ecal_file",type=str, required=True)
argparser.add_argument("--configs", help="configs", type=str, required=True)
argparser.add_argument("--plot_file", help="plot_file",type=str)
argparser.add_argument("--hit_pars", help="hit_pars",type=str)
args = argparser.parse_args()

with open(args.files[0]) as f:
    files = f.read().splitlines()
    files = sorted(files)

with open(args.ecal_file, 'r') as o:
    cal_dict = json.load(o)

energy_param = 'cuspEmax'
cal_energy_param = 'cuspEmax_ctc'
    

#out_dict = cal_aoe(files, cal_dict, energy_param, cal_energy_param, dt_corr=False, cut_parameters=cut_parameters, plot_savepath=args.plot_file)

out_dict = {"ecal_pars":[1,1],"aoe_pars":[1,1]}

if args.hit_pars is not None:
    pathlib.Path(os.path.dirname(args.hit_pars)).mkdir(parents=True, exist_ok=True)
    with open(args.hit_pars, 'w') as w:
        json.dump(out_dict,w, indent=4)
pathlib.Path(args.plot_file).touch()