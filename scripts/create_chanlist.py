# ruff: noqa: F821, T201

import os
import random
import re

from util.FileKey import ChannelProcKey
from util.patterns import get_pattern_pars_tmp_channel, get_pattern_plts_tmp_channel
from util.utils import filelist_path, runcmd

setup = snakemake.params.setup
basedir = snakemake.params.basedir
configs = snakemake.params.configs
chan_maps = snakemake.params.chan_maps

tier = snakemake.wildcards.tier
keypart = snakemake.wildcards.label

tier_pattern = (
    "((?P<file_type>[^_]+)(\\_(?P<tier>[^_]+)(\\_(?P<name>[^_]+)(\\_(?P<extension>[^_]+)?)?)?)?)?"
)
keypart_rx = re.compile(tier_pattern)
d = keypart_rx.match(tier).groupdict()

key = ChannelProcKey.parse_keypart(keypart)

output_file = os.path.join(
    filelist_path(setup),
    f"all-{key.experiment}-{key.period}-{key.run}-cal-{key.timestamp}-channels.chankeylist.{random.randint(0,99999):05d}",
)

cmd = f"{runcmd(setup)} python3 -B {basedir}/scripts/create_chankeylist.py --configs {configs}"
cmd += f" --channelmap {chan_maps} --timestamp {key.timestamp} --datatype cal --output_file {output_file}"
os.system(cmd)

with open(output_file) as r:
    chan_list = r.read().splitlines()

if d["file_type"] == "par":
    if d["extension"] is None:
        par_pattern = get_pattern_pars_tmp_channel(setup, d["tier"], d["name"])
    else:
        par_pattern = get_pattern_pars_tmp_channel(setup, d["tier"], d["name"], d["extension"])
elif d["file_type"] == "plt":
    par_pattern = get_pattern_plts_tmp_channel(setup, d["tier"], d["name"])
else:
    msg = "unknown file type: known are par and plt"
    raise ValueError(msg)

filenames = ChannelProcKey.get_channel_files(keypart, par_pattern, chan_list)

with open(snakemake.output[0], "w") as f:
    for fn in filenames:
        f.write(f"{fn}\n")

os.remove(output_file)
