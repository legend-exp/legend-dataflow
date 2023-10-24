import argparse
import logging
import os
import pathlib

import numpy as np
from legendmeta import LegendMetadata
from legendmeta.catalog import Props

argparser = argparse.ArgumentParser()
argparser.add_argument("input", help="input file", type=str)
argparser.add_argument("output", help="output file", type=str)
argparser.add_argument("blind_curve", help="blinding curves file", type=str)
argparser.add_argument("--datatype", help="Datatype", type=str, required=True)
argparser.add_argument("--timestamp", help="Timestamp", type=str, required=True)
argparser.add_argument("--configs", help="config file", type=str)
argparser.add_argument("--chan_maps", help="chan map", type=str)
argparser.add_argument("--log", help="log file", type=str)
args = argparser.parse_args()

os.makedirs(os.path.dirname(args.log), exist_ok=True)
logging.basicConfig(level=logging.INFO, filename=args.log, filemode="w")

pathlib.Path(os.path.dirname(args.output)).mkdir(parents=True, exist_ok=True)

blinding_curve = Props.read_from(args.blind_curve)

pathlib.Path(args.output).touch()