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

import logging
import time

from pygama.pargen.dsp_optimize import run_one_dsp
import pygama.pargen.energy_optimisation as om
import pygama.lgdo.lh5_store as lh5
import pygama.math.peak_fitting as pgf

def get_wf_indexes(sorted_indexs, n_events):
    out_list = []
    if isinstance(n_events, list):
        for i in range(len(n_events)):
            new_list = []
            for idx,entry in enumerate(sorted_indexs):
                if (entry >= np.sum(n_events[:i])) and (entry < np.sum(n_events[:i+1])):
                    new_list.append(idx)
            out_list.append(new_list)       
    else:
        for i in range(int(len(sorted_indexs)/n_events)):
            new_list = []
            for idx,entry in enumerate(sorted_indexs):
                if (entry >= i*n_events) and (entry < (i+1)*n_events):
                    new_list.append(idx)
            out_list.append(new_list)
    return out_list  

argparser = argparse.ArgumentParser()
#argparser.add_argument("files", help="files", nargs='*',type=str)
argparser.add_argument("--raw_filelist", help="raw_filelist", type=str)
argparser.add_argument("--decay_const", help="decay_const", type=str, required=True)
argparser.add_argument("--configs", help="configs", type=str, required=True)

argparser.add_argument("--log", help="log_file", type=str)

argparser.add_argument("--datatype", help="Datatype", type=str, required=True)
argparser.add_argument("--timestamp", help="Timestamp", type=str, required=True)
argparser.add_argument("--channel", help="Channel", type=str, required=True)

argparser.add_argument("--final_dsp_pars", help="final_dsp_pars", type=str, required=True)
argparser.add_argument("--alpha_dict", help="alpha_dict", type=str, required=True)
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


peaks_keV = np.array([583.191, 727.330, 860.564, 1620.5, 2614.553])
kev_widths = [(25,30), (25,35),(25,40),(25,55), (70,70)]
funcs = [pgf.gauss_step_pdf, pgf.radford_pdf, pgf.radford_pdf,pgf.radford_pdf,pgf.radford_pdf, pgf.radford_pdf]

final_idxs = []
kwarg_dicts_cusp = []
kwarg_dicts_trap = []
kwarg_dicts_zac = []
n_events=[]
for peak in peaks_keV:
    peak_idx = np.where(peaks_keV == peak)[0][0]
    func = funcs[peak_idx]
    kev_width = kev_widths[peak_idx]
    idxs = om.event_selection(raw_files, f'{args.channel}/raw', dsp_config, db_dict, peaks_keV, peak_idx, kev_width)
    final_idxs += idxs.tolist()
    n_events.append(len(idxs))
    kwarg_dicts_cusp.append({'parameter':'cuspEmax', 'func':pgf.extended_radford_pdf, 
                'gof_func':pgf.radford_pdf,'peak':peak, 'kev_width':kev_width})
    kwarg_dicts_zac.append({'parameter':'zacEmax', 'func':pgf.extended_radford_pdf, 
                'gof_func':pgf.radford_pdf,'peak':peak, 'kev_width':kev_width})
    kwarg_dicts_trap.append({'parameter':'trapEmax', 'func':pgf.extended_radford_pdf, 
                'gof_func':pgf.radford_pdf,'peak':peak, 'kev_width':kev_width})

sort_index = np.argsort(final_idxs)
idx_list = get_wf_indexes(sort_index, n_events)
idxs = np.array(sorted(final_idxs))
sto =lh5.LH5Store()

baseline = sto.read_object(f'{args.channel}/raw/baseline', raw_files,idx=idxs, n_rows = np.sum(n_events))[0]
waveforms = sto.read_object(f'{args.channel}/raw/waveform', raw_files,idx=idxs, n_rows = np.sum(n_events))[0]
tb_data = lh5.Table(col_dict = { 'waveform' : waveforms, 'baseline':baseline } )

t1 = time.time()
log.info(f'Data Loaded in {(t1-t0)/60} minutes')


kwarg_dict = [{'peak_dicts':kwarg_dicts_cusp, 'ctc_param':'QDrift', 'idx_list':idx_list, 'peaks_keV': peaks_keV},
              {'peak_dicts':kwarg_dicts_zac, 'ctc_param':'QDrift', 'idx_list':idx_list, 'peaks_keV': peaks_keV},
              {'peak_dicts':kwarg_dicts_trap, 'ctc_param':'QDrift', 'idx_list':idx_list, 'peaks_keV': peaks_keV}]


sample_x = np.array([[1,1],[1,5], [12,1], [12,5], [4,3], [8,3]])

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
    max_cusp =15
else:
    max_cusp = np.ceil(np.nanmax(sample_y_cusp)*2)
if np.isnan(sample_y_zac).all():
    max_trap =15
else:
    max_zac = np.ceil(np.nanmax(sample_y_zac)*2)
if np.isnan(sample_y_trap).all():
    max_trap =15
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
        

#n_events_cusp = [res["n_sig"] for res in results_cusp]
#allowed_n_sig_cusp = 


bopt_cusp = om.BayesianOptimizer(acq_func="lcb",batch_size=10)
bopt_cusp.add_dimension("cusp", "sigma", 1, 16, "us")
bopt_cusp.add_dimension("cusp", "flat", 1, 5, "us")

bopt_zac = om.BayesianOptimizer(acq_func="lcb",batch_size=10)
bopt_zac.add_dimension("zac", "sigma", 1, 16, "us")
bopt_zac.add_dimension("zac", "flat", 1, 5, "us")

bopt_trap = om.BayesianOptimizer(acq_func="lcb",batch_size=10)
bopt_trap.add_dimension("etrap", "rise", 1, 12, "us")
bopt_trap.add_dimension("etrap", "flat", 1, 5, "us")

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
                                                    nan_val = nan_vals, n_iter=20)

db_dict.update(out_param_dict)

t2 = time.time()

log.info(f'Optimiser finished in {(t2-t1)/60} minutes')

pathlib.Path(os.path.dirname(args.final_dsp_pars)).mkdir(parents=True, exist_ok=True)
with open(args.final_dsp_pars, 'w') as w:
    json.dump(db_dict, w, indent=4)

out_alpha_dict = {}
out_alpha_dict["cuspEmax_ctc"] = {"expression": "cuspEmax*(1+dt_eff*@a)",
                              "parameters":{"a":bopt_cusp.optimal_results["alpha"]}}

out_alpha_dict["zacEmax_ctc"] = {"expression": "zacEmax*(1+dt_eff*@a)",
                              "parameters":{"a":bopt_zac.optimal_results["alpha"]}}
    
out_alpha_dict["trapEmax_ctc"] = {"expression": "trapEmax*(1+dt_eff*@a)",
                              "parameters":{"a":bopt_trap.optimal_results["alpha"]}}


pathlib.Path(os.path.dirname(args.alpha_dict)).mkdir(parents=True, exist_ok=True)
with open(args.alpha_dict, 'w') as w:
    json.dump(out_alpha_dict, w, indent=4)


pathlib.Path(os.path.dirname(args.qbb_grid_path)).mkdir(parents=True, exist_ok=True)
with open(args.qbb_grid_path,"wb") as f:
    pkl.dump(optimisers, f)



def match_config(parameters, opt_dicts):
    """
    Matches config to parameters
    """
    out_dict = {}
    for opt_dict in opt_dicts:
        key = list(opt_dict.keys())[0]
        if key =='cusp':
            out_dict['cuspEmax'] = opt_dict
        elif key =='zac':
            out_dict['zacEmax'] = opt_dict
        elif key =='etrap':
            out_dict['trapEmax'] = opt_dict
    return out_dict

def load_all_grids(files, parameters):
    """
    #Loads in optimizer grids
    """
    grid_dict = {}
    for param in parameters:
        peak_grids = []
        for file in files:
            with open(file,"rb") as d:
                grid = pkl.load(d)
            peak_grids.append(grid[param])
        grid_dict[param] = peak_grids
    return grid_dict
"""

argparser = argparse.ArgumentParser()
argparser.add_argument("files", help="files", nargs='*',type=str)
argparser.add_argument("--raw_filelist", help="raw_filelist", type=str)
argparser.add_argument("--decay_const", help="decay_const", type=str, required=True)
argparser.add_argument("--configs", help="configs", type=str, required=True)

argparser.add_argument("--datatype", help="Datatype", type=str, required=True)
argparser.add_argument("--timestamp", help="Timestamp", type=str, required=True)
argparser.add_argument("--channel", help="Channel", type=str, required=True)

argparser.add_argument("--final_dsp_pars", help="final_dsp_pars", type=str, required=True)
argparser.add_argument("--qbb_grid_path", help="qbb_grid_path", type=str)

argparser.add_argument("--plot_save_path", help="plot_save_path", type=str, required=False)
args = argparser.parse_args()
    
cfg_file = os.path.join(args.configs, 'key_resolve.jsonl')
configs = config_catalog.get_config(cfg_file, args.configs, args.timestamp, args.datatype)
config_dict = configs['snakemake_rules']['pars_dsp_eopt']["inputs"]['processing_chain'][args.channel]
opt_json = configs['snakemake_rules']['pars_dsp_eopt']["inputs"]["optimiser_config"][args.channel]

with open(opt_json, 'r') as r:
    opt_dict = json.load(r)

peak_energies = np.array([])
for f in args.files:
    filename = os.path.basename(f)
    peak=filename.split('.p')[0].split("-")[6]
    peak_energy = float(peak)
    peak_energies = np.append(peak_energies, peak_energy)

parameters = ['cuspEmax'] #, 'zacEmax', 'trapEmax'
matched_configs = match_config(parameters, opt_dict)
peak_grids = load_all_grids(args.files, parameters)

if args.plot_save_path:
    energy_db_dict, fwhm_dict, qbb_grid = om.get_filter_params(peak_grids, matched_configs, 
                                                            peak_energies, parameters, save_path = args.plot_save_path)
else:
    energy_db_dict, fwhm_dict, qbb_grid = om.get_filter_params(peak_grids, matched_configs, 
                                                            peak_energies, parameters)



# FTP optimisation

with open(args.raw_filelist) as f:
    file_list = f.read().splitlines()

raw_file = sorted(file_list)

with open(args.decay_const, 'r') as t:
    db_dict = json.load(t)[args.channel]

db_dict.update(energy_db_dict)

wf_idxs = om.event_selection(raw_file, f'{args.channel}/raw', config_dict, db_dict, np.array([238.632,   583.191, 727.330, 860.564, 1620.5, 2614.553]), 5, (70,70))
parameters=['cuspEftp']

opt_json = configs['snakemake_rules']['pars_dsp_eopt']["inputs"]["optimiser_config_ftp"][args.channel]
print('Loaded configs')

with open(opt_json, 'r') as r:
    opt_dict = json.load(r)

out_grids = om.run_optimisation_multiprocessed(raw_file, opt_dict, config_dict, lh5_path = f'{args.channel}/raw',
            cuts = wf_idxs, db_dict = db_dict, 
            fom = om.fom_FWHM_fit,  n_events=10000,
            processes=30, parameter=parameters, func=pgf.extended_radford_pdf, gof_func=pgf.radford_pdf,
            peak=2614.553, kev_width=(70,70)) #kev_width=(50,50) pgf.radford_peak

ftp_dict = {}
for i in range(len(out_grids)):
    values = np.array([out_grids[i][j]['fwhm_o_max'] for j in range(len(out_grids[i]))])
    out_index = np.nanargmin(values)
    db_dict_entry1 = list(opt_dict[i].keys())[0]
    db_dict_entry2 = list(opt_dict[i][db_dict_entry1].keys())[0]
    opt_d = opt_dict[i][db_dict_entry1][db_dict_entry2]
    opt_vals = np.arange(opt_d['start'], opt_d['end'], opt_d['spacing'])
    #try:
    #    opt_vals = [ f'{val:.2f}*{opt_d["unit"]}' for val in opt_vals]
    #except:
    opt_vals = [ f'{val:.2f}' for val in opt_vals]
    ftp_dict[db_dict_entry1] = {db_dict_entry2:opt_vals[out_index]}

for key in list(db_dict.keys()):
    try:
        db_dict[key].update(ftp_dict[key]) 
    except:
        pass


#Save db dict values
pathlib.Path(os.path.dirname(args.final_dsp_pars)).mkdir(parents=True, exist_ok=True)
with open(args.final_dsp_pars, 'w') as w:
    json.dump(db_dict, w, indent=4)

pathlib.Path(os.path.dirname(args.qbb_grid_path)).mkdir(parents=True, exist_ok=True)
with open(args.qbb_grid_path,"wb") as f:
    pkl.dump({},f)
"""