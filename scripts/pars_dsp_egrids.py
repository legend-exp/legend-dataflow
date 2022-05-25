import json, os
import pygama.pargen.energy_optimising as om
from utils import run_splitter
import pygama.analysis.peak_fitting as pgf
from collections import OrderedDict
import pickle
import argparse
import pathlib
import time
import numpy as np

argparser = argparse.ArgumentParser()
argparser.add_argument("raw_files", help="raw_files", nargs='*',type=str)
argparser.add_argument("--output_path", help="output_path", type=str)
argparser.add_argument("--configs", help="configs", type=str, required=True)
argparser.add_argument("--decay_const", help="decay_const", type=str, required=True)
args = argparser.parse_args()

save_dict = {}
    
pathlib.Path(os.path.dirname(args.output_path)).mkdir(parents=True, exist_ok=True)
with open(args.output_path,"wb") as f:
    pickle.dump(save_dict,f)
