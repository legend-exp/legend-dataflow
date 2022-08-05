import snakemake as smk
import os, re, glob
import pygama.lgdo.lh5_store as lh5
import argparse
from util.metadata_loading import *
#from util.utils import *

#Placeholder chankeylist creator, will take in the tcm files and output the channels present


argparser = argparse.ArgumentParser()
argparser.add_argument("input", help="input files", nargs='*',type=str)

argparser.add_argument("--configs", help="configs", type=str, required=True)
argparser.add_argument("--datatype", help="Datatype", type=str, required=True)
argparser.add_argument("--timestamp", help="Timestamp", type=str, required=True)

argparser.add_argument("--output_file", help="output_file", type=str, required=True)
args = argparser.parse_args()

cfg_file = os.path.join(args.configs, 'key_resolve.jsonl')
configs = config_catalog.get_config(cfg_file, args.configs, args.timestamp, args.datatype)
channel_map = configs["hardware_configuration"]["channel_map"]

files = args.input
if isinstance(files, str):
    files= [files]

#channels = lh5.ls(files[0])
#for f in files[1:]:
#    chans = lh5.ls(f)
#    for chan in chans:
#        if chan not in channels:
#            channels.append(chan)

exp, per, run, datattype, tstamp, proc_step = os.path.basename(files[0]).split("-")

if run =="r006":

    channels = ["ch002", "ch003", "ch004", "ch005", "ch006", "ch013",
            "ch014",  "ch015",  "ch022", "ch023", "ch029", "ch030",
            "ch035",  "ch037", "ch038", "ch039", "ch040", "ch042", "ch043"]

elif run == "r012":

    channels = ["ch002", "ch003", "ch005", "ch006", "ch008", "ch009","ch010", "ch011", 
                "ch012", "ch013", "ch014", "ch021", "ch022", "ch023", "ch024", "ch025", "ch026",
                "ch027", "ch035", "ch036", "ch037", "ch038", "ch039", "ch040"]

elif run =="r014":

    channels = ["ch002", "ch003", "ch005", "ch006", "ch008", "ch009"," ch010", "ch011", 
                "ch013", "ch014", "ch016", "ch023", "ch024", "ch025", "ch026",
                "ch027", "ch028", "ch029", "ch037", "ch038", "ch039", "ch040", "ch041", "ch042", "ch043"]


#"ch001", "ch016", "ch036",

with open(args.output_file, 'w') as f:
    for chan in channels:
            f.write(f"{chan}\n")