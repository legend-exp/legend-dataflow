from __future__ import annotations

import argparse
import pickle as pkl
import shelve
from pathlib import Path

from dbetto.catalog import Props
from lgdo import lh5

from legenddataflow.methods import ChannelProcKey


def replace_path(d, old_path, new_path):
    if isinstance(d, dict):
        for k, v in d.items():
            d[k] = replace_path(v, old_path, new_path)
    elif isinstance(d, list):
        for i in range(len(d)):
            d[i] = replace_path(d[i], old_path, new_path)
    elif isinstance(d, str) and old_path in d:
        d = d.replace(old_path, new_path)
        d = d.replace(new_path, f"$_/{Path(new_path).name}")
    return d


def merge_channels() -> None:
    argparser = argparse.ArgumentParser()
    argparser.add_argument(
        "--input", help="input file", nargs="*", type=str, required=True
    )
    argparser.add_argument("--output", help="output file", type=str, required=True)
    argparser.add_argument(
        "--in-db",
        help="in db file (used for when lh5 files referred to in db)",
        type=str,
        required=False,
    )
    argparser.add_argument(
        "--out-db",
        help="lh5 file (used for when lh5 files referred to in db)",
        type=str,
        required=False,
    )
    args = argparser.parse_args()

    # change to only have 1 output file for multiple inputs
    # don't care about processing step, check if extension matches

    channel_files = args.input.infiles if hasattr(args.input, "infiles") else args.input

    file_extension = Path(args.output).suffix

    if file_extension in (".dat", ".dir"):
        common_dict = {}
        out_file = Path(args.output).with_suffix("")
    else:
        out_file = args.output

    Path(args.output).parent.mkdir(parents=True, exist_ok=True)

    if file_extension in (".json", ".yaml", ".yml"):
        out_dict = {}
        for channel in channel_files:
            if Path(channel).suffix == file_extension:
                channel_dict = Props.read_from(channel)
                fkey = ChannelProcKey.get_filekey_from_pattern(Path(channel).name)
                out_dict[fkey.channel] = channel_dict
            else:
                msg = "Output file extension does not match input file extension"
                raise RuntimeError(msg)

        Props.write_to(out_file, out_dict)

    elif file_extension == ".pkl":
        out_dict = {}
        for channel in channel_files:
            with Path(channel).open("rb") as r:
                channel_dict = pkl.load(r)
            fkey = ChannelProcKey.get_filekey_from_pattern(Path(channel).name)
            out_dict[fkey.channel] = channel_dict

        with Path(out_file).open("wb") as w:
            pkl.dump(out_dict, w, protocol=pkl.HIGHEST_PROTOCOL)

        Path(out_file).rename(out_file)

    elif file_extension in (".dat", ".dir"):
        common_dict = {}
        with shelve.open(str(out_file), "c", protocol=pkl.HIGHEST_PROTOCOL) as shelf:
            for channel in channel_files:
                with Path(channel).open("rb") as r:
                    channel_dict = pkl.load(r)
                fkey = ChannelProcKey.get_filekey_from_pattern(Path(channel).name)
                if isinstance(channel_dict, dict) and "common" in list(channel_dict):
                    chan_common_dict = channel_dict.pop("common")
                    common_dict[fkey.channel] = chan_common_dict
                shelf[fkey.channel] = channel_dict
            if len(common_dict) > 0:
                shelf["common"] = common_dict

    elif file_extension == ".lh5":
        if args.in_db:
            db_dict = Props.read_from(args.in_db)
        for channel in channel_files:
            if Path(channel).suffix == file_extension:
                fkey = ChannelProcKey.get_filekey_from_pattern(Path(channel).name)
                tb_in = lh5.read(f"{fkey.channel}", channel)

                lh5.write(
                    tb_in,
                    name=fkey.channel,
                    lh5_file=out_file,
                    wo_mode="a",
                )
                if args.in_db:
                    db_dict[fkey.channel] = replace_path(
                        db_dict[fkey.channel], channel, args.output
                    )
            else:
                msg = "Output file extension does not match input file extension"
                raise RuntimeError(msg)
        if args.out_db:
            Props.write_to(args.out_db, db_dict)
