import argparse, os, pathlib
import logging

import pygama
from pygama.raw.build_raw import * 
import numpy as np

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
  },
  "OrcaHeaderDecoder" : {
    "OrcaHeader" : {
      "key_list" : [ "*" ],
      "out_stream" : "{filekey}"
    }
  },
  "ORFlashCamListenerConfigDecoder" : {
    "FCConfig" : {
      "key_list" : [ "*" ],
      "out_stream" : "{filekey}"
    }
  },
  "*" : {
    "{name}" : {
      "key_list" : [ "*" ],
      "out_stream" : "{filekey}"
    }
  }
}
'''

rand_num = f'{np.random.randint(0,99999):05d}'
temp_output = f'{args.output}.{rand_num}'

build_raw(args.input, in_stream_type='ORCA', out_spec=json.loads(out_spec), filekey = temp_output, buffer_size=1024)

os.rename(temp_output, args.output)
