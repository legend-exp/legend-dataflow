import snakemake as smk
import os, re, glob
import pygama.lgdo.lh5_store as lh5
import argparse
from util.metadata_loading import *
#from util.utils import *

#Placeholder chankeylist creator, will take in the tcm files and output the channels present


argparser = argparse.ArgumentParser()
#argparser.add_argument("input", help="input files", nargs='*',type=str)

argparser.add_argument("--configs", help="configs", type=str, required=True)
argparser.add_argument("--datatype", help="Datatype", type=str, required=True)
argparser.add_argument("--timestamp", help="Timestamp", type=str, required=True)

argparser.add_argument("--output_file", help="output_file", type=str, required=True)
args = argparser.parse_args()

cfg_file = os.path.join(args.configs, 'key_resolve.jsonl')
configs = config_catalog.get_config(cfg_file, args.configs, args.timestamp, args.datatype)
channel_map = configs["hardware_configuration"]["channel_map"]

#files = args.input
#if isinstance(files, str):
#    files= [files]

#channels = lh5.ls(files[0])
#for f in files[1:]:
#    chans = lh5.ls(f)
#    for chan in chans:
#        if chan not in channels:
#            channels.append(chan)

channels = [chan for chan in channel_map if channel_map[chan]["software_status"]=="On" ]

with open(args.output_file, 'w') as f:
    for chan in channels:
            f.write(f"{chan}\n")