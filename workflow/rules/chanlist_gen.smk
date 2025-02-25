# ruff: noqa: F821, T201

import os
import random
import re

from legenddataflow.FileKey import ChannelProcKey
from legenddataflow.patterns import (
    get_pattern_pars_tmp_channel,
    get_pattern_plts_tmp_channel,
)
from legenddataflow import execenv_pyexe
from legenddataflow.utils import filelist_path


# FIXME: the system argument should always be explicitly supplied
def get_chanlist(
    setup, keypart, workflow, config, det_status, chan_maps, system="geds"
):
    key = ChannelProcKey.parse_keypart(keypart)

    flist_path = filelist_path(setup)
    os.makedirs(flist_path, exist_ok=True)
    output_file = os.path.join(
        flist_path,
        f"all-{key.experiment}-{key.period}-{key.run}-{key.datatype}-{key.timestamp}-channels.chankeylist.{random.randint(0,99999):05d}",
    )

    os.system(
        execenv_pyexe(config, "create-chankeylist")
        + f"--det-status {det_status} --channelmap {chan_maps} --timestamp {key.timestamp} "
        f"--datatype {key.datatype} --output-file {output_file} --system {system}"
    )

    with open(output_file) as r:
        chan_list = r.read().splitlines()
    os.remove(output_file)
    return chan_list


def get_par_chanlist(
    setup,
    keypart,
    tier,
    basedir,
    det_status,
    chan_maps,
    datatype="cal",
    system="geds",
    name=None,
    extension="yaml",
):

    chan_list = get_chanlist(
        setup, keypart, workflow, config, det_status, chan_maps, system
    )

    par_pattern = get_pattern_pars_tmp_channel(
        setup, tier, name, datatype=datatype, extension=extension
    )

    filenames = ChannelProcKey.get_channel_files(keypart, par_pattern, chan_list)

    return filenames


def get_plt_chanlist(setup, keypart, tier, basedir, det_status, chan_maps, name=None):

    chan_list = get_chanlist(setup, keypart, workflow, config, det_status, chan_maps)

    par_pattern = get_pattern_plts_tmp_channel(setup, tier, name)

    filenames = ChannelProcKey.get_channel_files(keypart, par_pattern, chan_list)

    return filenames
