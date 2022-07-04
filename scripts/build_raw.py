import argparse, os, pathlib
import logging

import pygama
from pygama.raw.build_raw import * 

argparser = argparse.ArgumentParser()
argparser.add_argument("input", help="input file", type=str)
argparser.add_argument("output", help="output file", type=str)
argparser.add_argument("--configs", help="config file", type=str)
argparser.add_argument("--log", help="log file", type=str)
args = argparser.parse_args()

logging.basicConfig(level=logging.INFO, filename=args.log, filemode='w')

pathlib.Path(os.path.dirname(args.output)).mkdir(parents=True, exist_ok=True)

# ToDo: Atomic file creation

out_spec = '''
{
  "ORFlashCamADCWaveformDecoder" : {
    "ch{key:03d}/raw" : {
      "key_list" : ["*"],
      "out_stream" : "{filekey}"
    }
  }
}
'''

build_raw(args.input, in_stream_type='ORCA', out_spec=json.loads(out_spec), filekey = args.output)
