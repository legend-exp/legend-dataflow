import json
from pathlib import Path

from dbetto import TextDB


def get_table_name(channelmap, timestamp, datatype, detector):
    if isinstance(channelmap, (str, Path)):
        channelmap = TextDB(channelmap, lazy=True)
    channel_dict = channelmap.valid_for(timestamp, system=datatype)
    return f"ch{channel_dict[detector].daq.rawid:07}"


def get_all_channels(channelmap, timestamp, datatype):
    if isinstance(channelmap, (str, Path)):
        channelmap = TextDB(channelmap, lazy=True)

    if isinstance(channelmap, TextDB):
        chmap = channelmap.channelmaps.on(timestamp, system=datatype)
    else:
        chmap = channelmap.valid_for(timestamp, system=datatype)

    channels = list(chmap)

    if len(channels) == 0:
        print("WARNING: No channels found")  # noqa: T201

    return channels


def get_table_mapping(channelmap, timestamp, datatype, tier):
    if isinstance(channelmap, (str, Path)):
        channelmap = TextDB(channelmap, lazy=True)
    channel_dict = channelmap.channelmap(timestamp, system=datatype)
    detectors = get_all_channels(channelmap, timestamp, datatype)
    return json.dumps(
        {
            detector: f"ch{channel_dict[detector].daq.rawid:07}/{tier}"
            for detector in detectors
        }
    )
