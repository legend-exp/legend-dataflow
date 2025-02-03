import argparse
from pathlib import Path

from dbetto import TextDB
from legendmeta import LegendMetadata


def create_chankeylist() -> None:
    argparser = argparse.ArgumentParser()
    argparser.add_argument("--det_status", help="det_status", type=str, required=True)
    argparser.add_argument("--datatype", help="Datatype", type=str, required=True)
    argparser.add_argument("--timestamp", help="Timestamp", type=str, required=True)
    argparser.add_argument("--channelmap", help="Channel Map", type=str, required=True)

    argparser.add_argument("--output_file", help="output_file", type=str, required=True)
    args = argparser.parse_args()

    det_status = TextDB(args.det_status, lazy=True)
    status_map = det_status.statuses.on(args.timestamp, system=args.datatype)

    channel_map = LegendMetadata(args.channelmap, lazy=True)
    chmap = channel_map.channelmaps.on(args.timestamp)

    channels = [
        chan
        for chan in status_map
        if status_map[chan]["processable"] is True and chmap[chan].system == "geds"
    ]
    Path(args.output_file).parent.mkdir(parents=True, exist_ok=True)
    with Path(args.output_file).open("w") as f:
        for chan in channels:
            f.write(f"{chan}\n")
