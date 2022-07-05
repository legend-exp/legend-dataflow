import json, os
import pygama.pargen.energy_optimisation as om
from util.utils import run_splitter
import pygama.math.peak_fitting as pgf
from collections import OrderedDict
import pickle
import argparse
import pathlib
import time
import numpy as np
from util.metadata_loading import *

argparser = argparse.ArgumentParser()
#argparser.add_argument("raw_files", help="raw_files", nargs='*',type=str)
argparser.add_argument("raw_files", help="raw_filelist", type=str)
argparser.add_argument("--configs", help="configs", type=str, required=True)
argparser.add_argument("--decay_const", help="decay_const", type=str, required=True)

argparser.add_argument("--datatype", help="Datatype", type=str, required=True)
argparser.add_argument("--timestamp", help="Timestamp", type=str, required=True)
argparser.add_argument("--channel", help="Channel", type=str, required=True)
argparser.add_argument("--peak", help="peak", type=float, required=True)

argparser.add_argument("--output_path", help="output_path", type=str)
args = argparser.parse_args()

peaks_keV = np.array([238.632,   583.191, 727.330, 860.564, 1620.5, 2614.553])
#kev_widths = [(10,10), (25,40), (25,40),(25,40),(25,40), (50,50)]
#funcs = [pgf.gauss_step, pgf.radford_peak, pgf.radford_peak,pgf.radford_peak,pgf.radford_peak, pgf.radford_peak]

if args.peak == 2614.553:
        kev_widths = (70, 70)
        n_processes = 19
        func = pgf.extended_radford_pdf
        gof_func = pgf.radford_pdf
elif args.peak == 238.632:
    kev_widths = (10,10)
    n_processes = 5
    func = pgf.extended_gauss_step_pdf
    gof_func = pgf.gauss_step_pdf
else:
    kev_widths = (25,55)
    n_processes = 19
    func = pgf.extended_radford_pdf
    gof_func = pgf.radford_pdf
    
peak_idx = np.where(peaks_keV == args.peak)[0][0]

with open(args.raw_files) as f:
    files = f.read().splitlines()

raw_files = sorted(files)

cfg_file = os.path.join(args.configs, 'key_resolve.jsonl')
configs = config_catalog.get_config(cfg_file, args.configs, args.timestamp, args.datatype)
config_dict = configs['snakemake_rules']['pars_dsp_egrid']["inputs"]['processing_chain'][args.channel]
opt_json = configs['snakemake_rules']['pars_dsp_egrid']["inputs"]["optimiser_config"][args.channel]

with open(args.decay_const, 'r') as t:
    db_dict = json.load(t)[args.channel]

with open(opt_json, 'r') as r:
    opt_dict = json.load(r)

wf_idxs = om.event_selection(raw_files, f'{args.channel}/raw', config_dict, db_dict, peaks_keV, peak_idx, kev_widths)


print('Loaded configs')

parameters=['cuspEmax'] #'zacEmax', 'trapEmax', 

t0 = time.time()

grid_out = om.run_optimisation_multiprocessed(raw_files, opt_dict, config_dict, cuts = wf_idxs, 
            lh5_path = f'{args.channel}/raw' , fom = om.fom_all_fit, db_dict = db_dict,  n_events=10000, 
            processes=n_processes, parameter=parameters, func=func, gof_func = gof_func,
            peak=args.peak, kev_width=kev_widths) 
    
t1 = time.time()

print(f'Calculated Grid in {(t1-t0)/60} minutes')

save_dict = {}
for i,param in enumerate(parameters):
    save_dict[param] = grid_out[i]
    
pathlib.Path(os.path.dirname(args.output_path)).mkdir(parents=True, exist_ok=True)
with open(args.output_path,"wb") as f:
    pickle.dump(save_dict,f)
