import json, os
import pygama.pargen.energy_optimisation as om
import pygama.math.peak_fitting as pgf
from collections import OrderedDict
import argparse
import pathlib
import pickle as pkl
import numpy as np
from util.metadata_loading import *

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
    Loads in optimizer grids
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
