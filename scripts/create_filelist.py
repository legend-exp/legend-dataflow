import snakemake as smk
import os, re, glob

from utils import *

setup = snakemake.params.setup

dataset_file = snakemake.input[0]
tier = snakemake.wildcards.tier

filenames = tier_files(setup, dataset_file, tier)
with open(snakemake.output[0], 'w') as f:
    for fn in filenames:
        f.write(f"{fn}\n")

