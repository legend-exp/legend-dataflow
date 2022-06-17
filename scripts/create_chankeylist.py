import snakemake as smk
import os, re, glob

#from util.utils import *

#Placeholder chankeylist creator, will take in the tcm files and output the channels present

channels = ["001", "002"]

with open(snakemake.output[0], 'w') as f:
    for chan in channels:
            f.write(f"{chan}\n")