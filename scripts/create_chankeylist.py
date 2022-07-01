import snakemake as smk
import os, re, glob
import pygama.lgdo.lh5_store as lh5
import argparse
#from util.utils import *

#Placeholder chankeylist creator, will take in the tcm files and output the channels present


argparser = argparse.ArgumentParser()
argparser.add_argument("--output_file", help="output_file", type=str, required=True)
argparser.add_argument("input", help="input files", nargs='*',type=str)
args = argparser.parse_args()

#files = args.input
#if isinstance(files, str):
#    files= [files]

#channels = lh5.ls(files[0])
#for f in files[1:]:
#    chans = lh5.ls(f)
#    for chan in chans:
#        if chan not in channels:
#            channels.append(chan)

channels = ["ch001", "ch002", "ch003", "ch004", "ch005", "ch006", "ch013",
           "ch014",  "ch015", "ch016", "ch022", "ch023", "ch029", "ch030",
           "ch035", "ch036", "ch037", "ch038", "ch039", "ch040", "ch042", "ch043"]


with open(args.output_file, 'w') as f:
    for chan in channels:
            f.write(f"{chan}\n")