
# ruff: noqa: F821, T201

import os
import random
import re

from scripts.util.FileKey import ChannelProcKey
from scripts.util.patterns import get_pattern_pars_tmp_channel, get_pattern_plts_tmp_channel
from scripts.util.utils import filelist_path, runcmd


def get_par_chanlist(setup, keypart, tier, basedir, configs, chan_maps, name=None, extension="json"):
    tier_pattern = (
        "((?P<file_type>[^_]+)(\\_(?P<tier>[^_]+)(\\_(?P<name>[^_]+)?)?)?)?"
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

    par_pattern = get_pattern_pars_tmp_channel(setup, tier, name, extension)

    filenames = ChannelProcKey.get_channel_files(keypart, par_pattern, chan_list)
    os.remove(output_file)
    return filenames

def get_plt_chanlist(setup, keypart, tier, basedir, configs, chan_maps, name=None):
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

    par_pattern = get_pattern_plts_tmp_channel(setup, tier, name)

    filenames = ChannelProcKey.get_channel_files(keypart, par_pattern, chan_list)
    os.remove(output_file)
    return filenames