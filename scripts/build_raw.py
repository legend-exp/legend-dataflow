import argparse, os, pathlib

import pygama
from pygama.io.daq_to_raw import * 

argparser = argparse.ArgumentParser()
argparser.add_argument("input", help="input file", type=str)
argparser.add_argument("output", help="output file", type=str)
argparser.add_argument("--configs", help="config file", type=str)
args = argparser.parse_args()

pathlib.Path(os.path.dirname(args.output)).mkdir(parents=True, exist_ok=True)

# ToDo: Atomic file creation

daq_to_raw(args.input, args.output, config={"daq":"ORCA"})
