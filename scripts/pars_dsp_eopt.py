from pygama.dsp.utils import numba_defaults

numba_defaults.cache = False
numba_defaults.boundscheck = True

import json, os
from collections import OrderedDict
import argparse
import pathlib
import pickle as pkl
import numpy as np
from util.metadata_loading import *
from sklearn.gaussian_process.kernels import *

import logging
import time

from pygama.pargen.dsp_optimize import run_one_dsp
import pygama.pargen.energy_optimisation as om
import pygama.lgdo.lh5_store as lh5
import pygama.math.peak_fitting as pgf

argparser = argparse.ArgumentParser()
argparser.add_argument("--raw_filelist", help="raw_filelist", type=str)
argparser.add_argument("--decay_const", help="decay_const", type=str, required=True)
argparser.add_argument("--configs", help="configs", type=str, required=True)

argparser.add_argument("--log", help="log_file", type=str)

argparser.add_argument("--datatype", help="Datatype", type=str, required=True)
argparser.add_argument("--timestamp", help="Timestamp", type=str, required=True)
argparser.add_argument("--channel", help="Channel", type=str, required=True)

argparser.add_argument("--final_dsp_pars", help="final_dsp_pars", type=str, required=True)
argparser.add_argument("--qbb_grid_path", help="qbb_grid_path", type=str)

argparser.add_argument("--plot_save_path", help="plot_save_path", type=str, required=False)
args = argparser.parse_args()
    
logging.basicConfig(level=logging.DEBUG, filename=args.log, filemode='w')
logging.getLogger('numba').setLevel(logging.INFO)
logging.getLogger('parse').setLevel(logging.INFO)
logging.getLogger('pygama.lgdo.lh5_store').setLevel(logging.INFO)
logging.getLogger('h5py._conv').setLevel(logging.INFO)
logging.getLogger('pygama.dsp.processing_chain').setLevel(logging.INFO)


log = logging.getLogger(__name__)


t0 = time.time()

cfg_file = os.path.join(args.configs, 'key_resolve.jsonl')
configs = config_catalog.get_config(cfg_file, args.configs, args.timestamp, args.datatype)
dsp_config = configs['snakemake_rules']['pars_dsp_eopt']["inputs"]['processing_chain'][args.channel]
opt_json = configs['snakemake_rules']['pars_dsp_eopt']["inputs"]["optimiser_config"][args.channel]

with open(opt_json, 'r') as r:
    opt_dict = json.load(r)

with open(args.raw_filelist) as f:
    files = f.read().splitlines()

raw_files = sorted(files)

with open(args.decay_const, 'r') as t:
    db_dict = json.load(t)

peaks_keV = np.array(opt_dict["peaks"])
kev_widths = [tuple(kev_width) for kev_width in opt_dict["kev_widths"]]


kwarg_dicts_cusp = []
kwarg_dicts_trap = []
kwarg_dicts_zac = []
for peak in peaks_keV:
    peak_idx = np.where(peaks_keV == peak)[0][0]
    kev_width = kev_widths[peak_idx]
    if peak == 238.632:
        kwarg_dicts_cusp.append({'parameter':'cuspEmax', 'func':pgf.extended_gauss_step_pdf, 
                    'gof_func':pgf.gauss_step_pdf,'peak':peak, 'kev_width':kev_width})
        kwarg_dicts_zac.append({'parameter':'zacEmax', 'func':pgf.extended_gauss_step_pdf, 
                    'gof_func':pgf.gauss_step_pdf,'peak':peak, 'kev_width':kev_width})
        kwarg_dicts_trap.append({'parameter':'trapEmax', 'func':pgf.extended_gauss_step_pdf, 
                    'gof_func':pgf.gauss_step_pdf,'peak':peak, 'kev_width':kev_width})
    else:
        
        kwarg_dicts_cusp.append({'parameter':'cuspEmax', 'func':pgf.extended_radford_pdf, 
                    'gof_func':pgf.radford_pdf,'peak':peak, 'kev_width':kev_width})
        kwarg_dicts_zac.append({'parameter':'zacEmax', 'func':pgf.extended_radford_pdf, 
                    'gof_func':pgf.radford_pdf,'peak':peak, 'kev_width':kev_width})
        kwarg_dicts_trap.append({'parameter':'trapEmax', 'func':pgf.extended_radford_pdf, 
                    'gof_func':pgf.radford_pdf,'peak':peak, 'kev_width':kev_width})

tb_data, idx_list = om.event_selection(raw_files, f'{args.channel}/raw', 
                                        dsp_config, db_dict,
                                         peaks_keV, np.arange(0,len(peaks_keV),1).tolist(),
                                         cut_parameters = opt_dict["cut_parameters"],
                                         n_events= opt_dict["n_events"]
                                         kev_widths)

t1 = time.time()
log.info(f'Data Loaded in {(t1-t0)/60} minutes')


kwarg_dict = [{'peak_dicts':kwarg_dicts_cusp, 'ctc_param':'QDrift', 'idx_list':idx_list, 'peaks_keV': peaks_keV},
              {'peak_dicts':kwarg_dicts_zac, 'ctc_param':'QDrift', 'idx_list':idx_list, 'peaks_keV': peaks_keV},
              {'peak_dicts':kwarg_dicts_trap, 'ctc_param':'QDrift', 'idx_list':idx_list, 'peaks_keV': peaks_keV}]


sample_x = np.array(opt_dict["initial_samples"])

results_cusp = []
results_zac = []
results_trap = []

sample_y_cusp = []
sample_y_zac = []
sample_y_trap = []

for i,x in enumerate(sample_x):
    
    db_dict["cusp"] = {"sigma":f'{x[0]}*us', "flat":f'{x[1]}*us'}
    db_dict["zac"] = {"sigma":f'{x[0]}*us', "flat":f'{x[1]}*us'}
    db_dict["etrap"] = {"rise":f'{x[0]}*us', "flat":f'{x[1]}*us'}

    log.info(f'Initialising values {i+1} : {db_dict}')
    
    tb_out = run_one_dsp(tb_data,
                dsp_config,
                db_dict=db_dict,
                verbosity=0)
    
    res = om.new_fom(tb_out, kwarg_dict[0])
    results_cusp.append(res)
    sample_y_cusp.append(res['y_val'])
    
    res = om.new_fom(tb_out, kwarg_dict[1])
    results_zac.append(res)
    sample_y_zac.append(res['y_val'])
    
    res = om.new_fom(tb_out, kwarg_dict[2])
    results_trap.append(res)
    sample_y_trap.append(res['y_val'])
    
    log.info(f'{i+1} Finished')

if np.isnan(sample_y_cusp).all():
    max_cusp =opt_dict["nan_default"]
else:
    max_cusp = np.ceil(np.nanmax(sample_y_cusp)*2)
if np.isnan(sample_y_zac).all():
    max_zac =opt_dict["nan_default"]
else:
    max_zac = np.ceil(np.nanmax(sample_y_zac)*2)
if np.isnan(sample_y_trap).all():
    max_trap =opt_dict["nan_default"]
else:
    max_trap = np.ceil(np.nanmax(sample_y_trap)*2)

nan_vals = [max_cusp, max_zac, max_trap]

for i in range(len(sample_x)):

    if np.isnan(sample_y_cusp[i]):
        results_cusp[i]['y_val']=max_cusp
        sample_y_cusp[i] = max_cusp
    
    if np.isnan(sample_y_zac[i]):
        results_zac[i]['y_val']=max_zac
        sample_y_zac[i] = max_zac
        
    if np.isnan(sample_y_trap[i]):
        results_trap[i]['y_val']=max_trap
        sample_y_trap[i] = max_trap


bopt_cusp = om.BayesianOptimizer(acq_func=opt_dict["acq_func"],batch_size=opt_dict["batch_size"])
bopt_cusp.kernel = ConstantKernel(1.9, constant_value_bounds="fixed")+ConstantKernel(1.0, constant_value_bounds="fixed") * RBF(1, length_scale_bounds="fixed")
bopt_cusp.add_dimension("cusp", "sigma", 1, 16, "us")
bopt_cusp.add_dimension("cusp", "flat", 1.5, 2.5, "us")

bopt_zac = om.BayesianOptimizer(acq_func=opt_dict["acq_func"],batch_size=opt_dict["batch_size"])
bopt_zac.kernel = ConstantKernel(1.9, constant_value_bounds="fixed")+ConstantKernel(1.0, constant_value_bounds="fixed") * RBF(1, length_scale_bounds="fixed")
bopt_zac.add_dimension("zac", "sigma", 1, 16, "us")
bopt_zac.add_dimension("zac", "flat", 1.5, 2.5, "us")

bopt_trap = om.BayesianOptimizer(acq_func=opt_dict["acq_func"],batch_size=opt_dict["batch_size"])
bopt_trap.kernel = ConstantKernel(1.9, constant_value_bounds="fixed")+ConstantKernel(1.0, constant_value_bounds="fixed") * RBF(1, length_scale_bounds="fixed")
bopt_trap.add_dimension("etrap", "rise", 1, 12, "us")
bopt_trap.add_dimension("etrap", "flat", 1.5, 2.5, "us")

bopt_cusp.add_initial_values(x_init=sample_x, y_init=sample_y_cusp)
bopt_zac.add_initial_values(x_init=sample_x, y_init=sample_y_zac)
bopt_trap.add_initial_values(x_init=sample_x, y_init=sample_y_trap)

best_idx = np.nanargmin(sample_y_cusp)
bopt_cusp.optimal_results = results_cusp[best_idx]
bopt_cusp.optimal_x = sample_x[best_idx]

best_idx = np.nanargmin(sample_y_zac)
bopt_zac.optimal_results = results_zac[best_idx]
bopt_zac.optimal_x = sample_x[best_idx]

best_idx = np.nanargmin(sample_y_trap)
bopt_trap.optimal_results = results_trap[best_idx]
bopt_trap.optimal_x = sample_x[best_idx]

optimisers = [bopt_cusp,bopt_zac,bopt_trap]

out_param_dict, out_results_list = om.run_optimisation(tb_data, dsp_config, [om.new_fom], optimisers, 
                                                    fom_kwargs=kwarg_dict,db_dict=db_dict, 
                                                    nan_val = nan_vals, n_iter=opt_dict["n_iter"])

db_dict.update(out_param_dict)

t2 = time.time()
log.info(f'Optimiser finished in {(t2-t1)/60} minutes')

out_alpha_dict = {}
out_alpha_dict["cuspEmax_ctc"] = {"expression": "cuspEmax*(1+dt_eff*@a)",
                              "parameters":{"a":bopt_cusp.optimal_results["alpha"]}}

out_alpha_dict["zacEmax_ctc"] = {"expression": "zacEmax*(1+dt_eff*@a)",
                              "parameters":{"a":bopt_zac.optimal_results["alpha"]}}
    
out_alpha_dict["trapEmax_ctc"] = {"expression": "trapEmax*(1+dt_eff*@a)",
                              "parameters":{"a":bopt_trap.optimal_results["alpha"]}}

db_dict.update({"ctc_params":out_alpha_dict})

pathlib.Path(os.path.dirname(args.final_dsp_pars)).mkdir(parents=True, exist_ok=True)
with open(args.final_dsp_pars, 'w') as w:
    json.dump(db_dict, w, indent=4)


pathlib.Path(os.path.dirname(args.qbb_grid_path)).mkdir(parents=True, exist_ok=True)
with open(args.qbb_grid_path,"wb") as f:
    pkl.dump(optimisers, f)
