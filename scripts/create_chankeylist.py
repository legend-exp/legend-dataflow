import snakemake as smk
import os, re, glob
import pygama.lgdo.lh5_store as lh5
import argparse
from legendmeta import LegendMetadata


argparser = argparse.ArgumentParser()
argparser.add_argument("--configs", help="configs", type=str, required=True)
argparser.add_argument("--datatype", help="Datatype", type=str, required=True)
argparser.add_argument("--timestamp", help="Timestamp", type=str, required=True)

argparser.add_argument("--output_file", help="output_file", type=str, required=True)
args = argparser.parse_args()

configs = LegendMetadata(path = args.configs)
channel_map = configs.on(args.timestamp, system=args.datatype)["hardware_configuration"]["channel_map"]


channels = [chan for chan in channel_map if channel_map[chan]["software_status"]=="On" ]

with open(args.output_file, 'w') as f:
    for chan in channels:
            f.write(f"{chan}\n")