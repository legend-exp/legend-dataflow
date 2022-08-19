import argparse, pathlib
import numpy as np
import os,json

from pygama.hit.build_hit import build_hit

import time
import logging

argparser = argparse.ArgumentParser()
argparser.add_argument("input", help="input file", type=str)
argparser.add_argument("--pars_file", help="hit pars file", type=str)

argparser.add_argument("--configs", help="configs", type=str, required=True)
argparser.add_argument("--datatype", help="Datatype", type=str, required=True)
argparser.add_argument("--timestamp", help="Timestamp", type=str, required=True)

argparser.add_argument("--log", help="log_file", type=str)

argparser.add_argument("--output", help="output file", type=str)
args = argparser.parse_args()

logging.basicConfig(level=logging.DEBUG, filename=args.log, filemode='w')
logging.getLogger('numba').setLevel(logging.INFO)
logging.getLogger('parse').setLevel(logging.INFO)
logging.getLogger('pygama.lgdo.lh5_store').setLevel(logging.INFO)
logging.getLogger('h5py._conv').setLevel(logging.INFO)


log = logging.getLogger(__name__)

with open(args.pars_file) as f:
    pars_dict = json.load(f)

hit_dict ={}
for channel in pars_dict:
    hit_dict[f'{channel}/dsp']=pars_dict[channel]

t_start = time.time()
pathlib.Path(os.path.dirname(args.output)).mkdir(parents=True, exist_ok=True)
build_hit(args.input, lh5_tables_config=hit_dict, outfile =args.output)
t_elap = time.time() - t_start
log.info(f'Done!  Time elapsed: {t_elap:.2f} sec.')