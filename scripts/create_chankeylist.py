import snakemake as smk
import os, re, glob

from utils import *

#Placeholder chankeylist creator, will take in the tcm files and output the channels present

channels = ["1", "2"]

with open(snakemake.output[0], 'w') as f:
    for chan in channels:
            f.write(f"{chan}\n")