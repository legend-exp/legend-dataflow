import json, os
import pygama.pargen.energy_optimising as om
from utils import run_splitter
import pygama.analysis.peak_fitting as pgf
from collections import OrderedDict
import argparse
import pathlib
import pickle as pkl
import numpy as np

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
argparser.add_argument("--final_dsp_pars", help="final_dsp_pars", type=str, required=True)
argparser.add_argument("--configs", help="configs", type=str, required=True)
argparser.add_argument("--qbb_grid_path", help="qbb_grid_path", type=str)
argparser.add_argument("--raw_filelist", help="raw_filelist", type=str)
argparser.add_argument("--decay_const", help="decay_const", type=str, required=True)
argparser.add_argument("--plot_save_path", help="plot_save_path", type=str, required=False)
args = argparser.parse_args()
    
final_db_dict = {
    "pz": {
	"tau": 2886.0
    },
    "zac": {
	"sigma": "7.0*us",
        "flat": "2.25*us",
        "alpha": 1.2070354256073482e-06,
        "sample": "40.00"
    },
    "etrap": {
	"rise": "7.0*us",
        "flat": "2.25*us",
        "alpha": 1.1846359723793476e-06,
        "sample": "0.80"
    },
    "cusp": {
	"sigma": "3.5*us",
        "flat": "2.25*us",
        "alpha": 1.1924692049451807e-06,
        "sample": "35.00"
    }
}

#Save db dict values
pathlib.Path(os.path.dirname(args.final_dsp_pars)).mkdir(parents=True, exist_ok=True)
with open(args.final_dsp_pars, 'w') as w:
    json.dump(final_db_dict, w, indent=4)

pathlib.Path(os.path.dirname(args.qbb_grid_path)).mkdir(parents=True, exist_ok=True)
with open(args.qbb_grid_path,"wb") as f:
    pkl.dump({},f)

pathlib.Path(args.plot_save_path).touch()