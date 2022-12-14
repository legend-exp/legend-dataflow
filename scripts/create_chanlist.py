import snakemake as smk
from util.FileKey import *
from util.patterns import *

setup = snakemake.params.setup

dataset_file = snakemake.input[0]
tier = snakemake.wildcards.tier
keypart = snakemake.wildcards.label
name=None
if "_" in tier:
    name = tier.split("_", 1)[1]
    tier = tier.split("_", 1)[0]
    

filenames = ChannelFileKey.get_channel_files(setup, keypart, tier, dataset_file, name=name)

with open(snakemake.output[0], 'w') as f:
    for fn in filenames:
        f.write(f"{fn}\n")

