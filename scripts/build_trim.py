import argparse, os, pathlib
import logging

import pygama
from data_trimmer.data_trimmer import data_trimmer
import numpy as np

argparser = argparse.ArgumentParser()
argparser.add_argument("input", help="input file", type=str)
argparser.add_argument("output", help="output file", type=str)
argparser.add_argument("--configs", help="config file", type=str)
argparser.add_argument("--log", help="log file", type=str)
args = argparser.parse_args()

logging.basicConfig(level=logging.INFO, filename=args.log, filemode='w')

# both trimmed files are sent to the same directory as of now
pathlib.Path(os.path.dirname(args.output[0])).mkdir(parents=True, exist_ok=True)

# TODO: move this to the actual config... and drop in the correct values for the processors
# Note: the window value is actually set in the data_trimmer.py script
trim_config = '''
{
    "outputs" : ["presummed" ],
    "processors" : {        
        "presummed": {
            "function": "presum",
            "module": "pygama.dsp.processors",
            "args": ["waveform", "presummed(len(waveform)/4, 'f')"],
            "unit": "ADC"
        }
    }
 }
'''

# the order of the arguments is: input file, windowed, presummed, dsp_config
data_trimmer(args.input, args.output[0], args.output[1], trim_config)
