import os
import random

from util.FileKey import *
from util.patterns import *

setup = snakemake.params.setup
basedir = snakemake.params.basedir
configs = snakemake.params.configs
chan_maps = snakemake.params.chan_maps

tier = snakemake.wildcards.tier
keypart = snakemake.wildcards.label
name = None
if "_" in tier:
    name = tier.split("_", 1)[1]
    tier = tier.split("_", 1)[0]

key = ChannelProcKey.parse_keypart(keypart)

output_file = os.path.join(
    filelist_path(setup),
    f"all-{key.experiment}-{key.period}-{key.run}-cal-{key.timestamp}-channels.chankeylist.{random.randint(0,99999):05d}",
)

cmd = f"{runcmd(setup)} python3 -B {basedir}/scripts/create_chankeylist.py --configs {configs} --channelmap {chan_maps} --timestamp {key.timestamp} --datatype cal --output_file {output_file}"
os.system(cmd)

with open(output_file) as r:
    chan_list = r.read().splitlines()

filenames = ChannelProcKey.get_channel_files(setup, keypart, tier, chan_list, name=name)

with open(snakemake.output[0], "w") as f:
    for fn in filenames:
        f.write(f"{fn}\n")

os.remove(output_file)
