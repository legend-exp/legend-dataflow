import snakemake as smk
import os, re, glob

from utils import *

setup = snakemake.params.setup

dataset_file = snakemake.input[0]
tier = snakemake.wildcards.tier
keypart = snakemake.wildcards.label

name = f'{tier}_pars'

d = parse_channel_keypart(keypart)

par_pattern = get_pattern_pars_tmp_channel(setup,tier,name)
par_pattern_rx = re.compile(smk.io.regex(par_pattern))

filenames = []
with open(dataset_file) as f:
    for channel in f:
        file = smk.io.expand(par_pattern, experiment = d["experiment"], 
                                period = d["period"], run = d["run"], datatype = d["datatype"], channel=channel.strip())[0]
        filenames.append(file)

with open(snakemake.output[0], 'w') as f:
    for fn in filenames:
        f.write(f"{fn}\n")