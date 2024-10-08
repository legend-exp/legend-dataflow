import argparse
import os
import pathlib

from legendmeta import LegendMetadata, TextDB

argparser = argparse.ArgumentParser()
argparser.add_argument("--configs", help="configs", type=str, required=True)
argparser.add_argument("--datatype", help="Datatype", type=str, required=True)
argparser.add_argument("--timestamp", help="Timestamp", type=str, required=True)
argparser.add_argument("--channelmap", help="Channel Map", type=str, required=True)

argparser.add_argument("--output_file", help="output_file", type=str, required=True)
args = argparser.parse_args()

configs = TextDB(args.configs, lazy=True)
status_map = configs.on(args.timestamp, system=args.datatype)["analysis"]

channel_map = LegendMetadata(args.channelmap, lazy=True)
chmap = channel_map.channelmaps.on(args.timestamp)

channels = [
    f"ch{chmap[chan].daq.rawid:03}"
    for chan in status_map
    if status_map[chan]["processable"] is True and chmap[chan].system == "geds"
]

pathlib.Path(os.path.dirname(args.output_file)).mkdir(parents=True, exist_ok=True)
with open(args.output_file, "w") as f:
    for chan in channels:
        f.write(f"{chan}\n")
