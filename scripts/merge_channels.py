import argparse
import os
import pathlib
import pickle as pkl
import shelve

import numpy as np
from legendmeta.catalog import Props
from lgdo import lh5
from util.FileKey import ChannelProcKey
from util.utils import as_ro


def replace_path(d, old_path, new_path):
    if isinstance(d, dict):
        for k, v in d.items():
            d[k] = replace_path(v, old_path, new_path)
    elif isinstance(d, list):
        for i in range(len(d)):
            d[i] = replace_path(d[i], old_path, new_path)
    elif isinstance(d, str) and old_path in d:
        d = d.replace(old_path, new_path)
        d = d.replace(new_path, f"$_/{os.path.basename(new_path)}")
    return d


argparser = argparse.ArgumentParser()
argparser.add_argument("--input", help="input file", nargs="*", type=str, required=True)
argparser.add_argument("--output", help="output file", type=str, required=True)
argparser.add_argument(
    "--in_db",
    help="in db file (used for when lh5 files referred to in db)",
    type=str,
    required=False,
)
argparser.add_argument(
    "--out_db",
    help="lh5 file (used for when lh5 files referred to in db)",
    type=str,
    required=False,
)
args = argparser.parse_args()

# change to only have 1 output file for multiple inputs
# don't care about processing step, check if extension matches
channel_files = as_ro(args.input)

file_extension = pathlib.Path(args.output).suffix

if file_extension == ".dat" or file_extension == ".dir":
    out_file = os.path.splitext(args.output)[0]
else:
    out_file = args.output

rng = np.random.default_rng()
temp_output = f"{out_file}.{rng.integers(0, 99999):05d}"

pathlib.Path(os.path.dirname(args.output)).mkdir(parents=True, exist_ok=True)

if file_extension == ".json" or file_extension == ".yaml" or file_extension == ".yml":
    out_dict = {}
    for channel in channel_files:
        if pathlib.Path(channel).suffix == file_extension:
            channel_dict = Props.read_from(channel)

            fkey = ChannelProcKey.get_filekey_from_pattern(os.path.basename(channel))
            channel_name = fkey.channel
            out_dict[channel_name] = channel_dict
        else:
            msg = "Output file extension does not match input file extension"
            raise RuntimeError(msg)

    Props.write_to(temp_output, out_dict)

    os.rename(temp_output, out_file)

elif file_extension == ".pkl":
    out_dict = {}
    for channel in channel_files:
        with open(channel, "rb") as r:
            channel_dict = pkl.load(r)
        fkey = ChannelProcKey.get_filekey_from_pattern(os.path.basename(channel))
        channel_name = fkey.channel
        out_dict[channel_name] = channel_dict

    with open(temp_output, "wb") as w:
        pkl.dump(out_dict, w, protocol=pkl.HIGHEST_PROTOCOL)

    os.rename(temp_output, out_file)

elif file_extension == ".dat" or file_extension == ".dir":
    common_dict = {}
    with shelve.open(out_file, "c", protocol=pkl.HIGHEST_PROTOCOL) as shelf:
        for channel in channel_files:
            with open(channel, "rb") as r:
                channel_dict = pkl.load(r)
            fkey = ChannelProcKey.get_filekey_from_pattern(os.path.basename(channel))
            channel_name = fkey.channel
            if isinstance(channel_dict, dict) and "common" in list(channel_dict):
                chan_common_dict = channel_dict.pop("common")
                common_dict[channel_name] = chan_common_dict
            shelf[channel_name] = channel_dict
        if len(common_dict) > 0:
            shelf["common"] = common_dict


elif file_extension == ".lh5":
    if args.in_db:
        db_dict = Props.read_from(args.in_db)
    for channel in channel_files:
        if pathlib.Path(channel).suffix == file_extension:
            fkey = ChannelProcKey.get_filekey_from_pattern(os.path.basename(channel))
            channel_name = fkey.channel

            tb_in = lh5.read(f"{channel_name}", channel)[0]

            lh5.write(
                tb_in,
                name=channel_name,
                lh5_file=temp_output,
                wo_mode="a",
            )
            if args.in_db:
                db_dict[channel_name] = replace_path(db_dict[channel_name], channel, args.output)
        else:
            msg = "Output file extension does not match input file extension"
            raise RuntimeError(msg)
    if args.out_db:
        Props.write_to(args.out_db, db_dict)

    os.rename(temp_output, out_file)
