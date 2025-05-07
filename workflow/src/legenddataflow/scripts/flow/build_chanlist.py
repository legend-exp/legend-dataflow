# ruff: noqa: F821, T201
from __future__ import annotations

from pathlib import Path

from dbetto import TextDB

from legenddataflow.methods import ChannelProcKey
from legenddataflow.methods.patterns import (
    get_pattern_pars_tmp_channel,
    get_pattern_plts_tmp_channel,
)


def get_chanlist(timestamp, datatype, det_status, channelmap, system):
    if isinstance(det_status, str | Path):
        det_status = TextDB(det_status, lazy=True)

    if isinstance(channelmap, str | Path):
        channelmap = TextDB(channelmap, lazy=True)

    if isinstance(det_status, TextDB):
        status_map = det_status.statuses.on(timestamp, system=datatype)
    else:
        status_map = det_status.valid_for(timestamp, system=datatype)
    if isinstance(channelmap, TextDB):
        chmap = channelmap.channelmaps.on(timestamp, system=datatype)
    else:
        chmap = channelmap.valid_for(timestamp, system=datatype)

    # only restrict to a certain system (geds, spms, ...)
    channels = []
    for channel in chmap.map("system", unique=False)[system].map("name"):
        if channel not in status_map:
            msg = f"{channel} is not found in the status map (on {timestamp})"
            raise RuntimeError(msg)
        if status_map[channel].processable is False:
            continue
        channels.append(channel)

    if len(channels) == 0:
        print("WARNING: No channels found")

    return channels


def get_par_chanlist(
    setup,
    keypart,
    tier,
    det_status,
    chan_maps,
    system,
    datatype="cal",
    name=None,
    extension="yaml",
):
    key = ChannelProcKey.parse_keypart(keypart)

    chan_list = get_chanlist(key.timestamp, key.datatype, det_status, chan_maps, system)

    par_pattern = get_pattern_pars_tmp_channel(
        setup, tier, name, datatype=datatype, extension=extension
    )

    return ChannelProcKey.get_channel_files(keypart, par_pattern, chan_list)


def get_plt_chanlist(
    setup,
    keypart,
    tier,
    det_status,
    chan_maps,
    system,
    name=None,
):
    key = ChannelProcKey.parse_keypart(keypart)

    chan_list = get_chanlist(key.timestamp, key.datatype, det_status, chan_maps, system)

    par_pattern = get_pattern_plts_tmp_channel(setup, tier, name)

    return ChannelProcKey.get_channel_files(keypart, par_pattern, chan_list)
