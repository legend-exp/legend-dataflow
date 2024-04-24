import argparse
import json
import os
import pathlib
import pickle as pkl
import shelve

import lgdo.lh5 as lh5
import numpy as np
from legendmeta.catalog import Props
from util.FileKey import ChannelProcKey


def replace_path(d, old_path, new_path):
    if isinstance(d, dict):
        for k, v in d.items():
            d[k] = replace_path(v, old_path, new_path)
    elif isinstance(d, list):
        for i in range(len(d)):
            d[i] = replace_path(d[i], old_path, new_path)
    elif isinstance(d, str) and old_path in d:
        d = d.replace(old_path, new_path)
        d = f"$_/{os.path.basename(new_path)}"
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


channel_files = args.input

file_extension = pathlib.Path(args.output).suffix

if file_extension == ".dat" or file_extension == ".dir":
    out_file = os.path.splitext(args.output)[0]
else:
    out_file = args.output

rng = np.random.default_rng()
rand_num = f"{rng.integers(0,99999):05d}"
temp_output = f"{out_file}.{rand_num}"

pathlib.Path(os.path.dirname(args.output)).mkdir(parents=True, exist_ok=True)


if file_extension == ".json":
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

    with open(temp_output, "w") as w:
        json.dump(out_dict, w, indent=4)

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
    sto = lh5.LH5Store()

    if args.in_db:
        db_dict = Props.read_from(args.in_db)
    for channel in channel_files:
        if pathlib.Path(channel).suffix == file_extension:
            fkey = ChannelProcKey.get_filekey_from_pattern(os.path.basename(channel))
            channel_name = fkey.channel

            tb_in = sto.read(f"{channel_name}", channel)[0]

            sto.write(
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
        with open(args.out_db, "w") as w:
            json.dump(db_dict, w, indent=4)

    os.rename(temp_output, out_file)
