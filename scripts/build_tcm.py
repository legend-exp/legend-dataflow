import  os, logging
import argparse, pathlib
from pygama.evt.build_tcm import *
from pygama import lgdo
from pygama.raw.orca import orca_flashcam
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

# get the list of channels by fcid
ch_list = lgdo.ls(args.input, '/ch*')
fcid_channels = {}
for ch in ch_list:
    key = int(ch[2:])
    fcid = orca_flashcam.get_fcid(key)
    if fcid not in fcid_channels: fcid_channels[fcid] = []
    fcid_channels[fcid].append(f'/{ch}/raw')

# make a hardware_tcm_[fcid] for each fcid
for fcid in fcid_channels.keys():
    out_name = f'hardware_tcm_{fcid}'
    ch_list = fcid_channels[fcid]
    build_tcm([(args.input, ch_list)], 'eventnumber', out_file=temp_output, out_name=out_name, wo_mode='o')

os.rename(temp_output, args.output)
