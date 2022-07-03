import  os, logging
import argparse, pathlib
from pygama.evt.build_tcm import *

argparser = argparse.ArgumentParser()
argparser.add_argument("input", help="input file", type=str)
argparser.add_argument("output", help="output file", type=str)
argparser.add_argument("--configs", help="config file", type=str)
argparser.add_argument("--log", help="log file", type=str)
args = argparser.parse_args()

logging.basicConfig(level=logging.DEBUG, filename=args.log, filemode='w')

print(args.input, args.output)

pathlib.Path(os.path.dirname(args.output)).mkdir(parents=True, exist_ok=True)

build_tcm([(args.input, '/ch*/raw')], 'eventnumber', out_file=args.output, out_name='hardware_tcm', wo_mode='o')