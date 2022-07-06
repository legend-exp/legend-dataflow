import  os, logging
import argparse, pathlib
from pygama.evt.build_tcm import *
import numpy as np

argparser = argparse.ArgumentParser()
argparser.add_argument("input", help="input file", type=str)
argparser.add_argument("output", help="output file", type=str)
argparser.add_argument("--configs", help="config file", type=str)
argparser.add_argument("--log", help="log file", type=str)
args = argparser.parse_args()

logging.basicConfig(level=logging.DEBUG, filename=args.log, filemode='w')

pathlib.Path(os.path.dirname(args.output)).mkdir(parents=True, exist_ok=True)

rand_num = f'{np.random.randint(0,99999):05d}'
temp_output = f'{args.output}.{rand_num}'

build_tcm([(args.input, '/ch*/raw')], 'eventnumber', out_file=temp_output, out_name='hardware_tcm', wo_mode='o')

os.rename(temp_output, args.output)