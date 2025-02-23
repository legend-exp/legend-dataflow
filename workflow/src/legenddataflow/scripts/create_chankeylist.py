import argparse
from pathlib import Path

from dbetto import TextDB


def create_chankeylist() -> None:
    argparser = argparse.ArgumentParser()

    argparser.add_argument("--det-status", help="det_status", required=True)
    argparser.add_argument("--datatype", help="Datatype", required=True)
    argparser.add_argument("--timestamp", help="Timestamp", required=True)
    argparser.add_argument("--channelmap", help="Channel Map", required=True)
    argparser.add_argument("--system", help="geds, spms, pmts, ...", required=True)
    argparser.add_argument("--output-file", required=True)

    args = argparser.parse_args()

    status_map = TextDB(args.det_status, lazy=True).statuses.on(
        args.timestamp, system=args.datatype
    )
    chmap = TextDB(args.channelmap, lazy=True).channelmaps.on(args.timestamp)

    # only restrict to a certain system (geds, spms, ...)
    channels = []
    for channel, status in status_map.items():
        # start with channels marked as processable in the status map
        if status.processable is False:
            continue

        if channel not in chmap:
            msg = f"{channel} is marked as processable but is not found in the channel map (on {args.timestamp})"
            raise RuntimeError(msg)

        if chmap[channel].system == args.system:
            channels.append(channel)

    if len(channels) == 0:
        print("WARNING: No channels found")  # noqa: T201

    Path(args.output_file).parent.mkdir(parents=True, exist_ok=True)
    with Path(args.output_file).open("w") as f:
        for chan in channels:
            f.write(f"{chan}\n")
