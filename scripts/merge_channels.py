import argparse
import json
import os
import pathlib
import pickle as pkl
import shelve

argparser = argparse.ArgumentParser()
argparser.add_argument("--input", help="input file", nargs="*", type=str)
argparser.add_argument("--output", help="output file", nargs="*", type=str)
args = argparser.parse_args()


channel_files = args.input
for _i, out_file in enumerate(args.output):
    file_extension = pathlib.Path(out_file).suffix
    processing_step = os.path.splitext(out_file)[0].split("-")[-1]
    if file_extension == ".json":
        out_dict = {}
        for channel in channel_files:
            if os.path.splitext(channel)[0].split("-")[-1] == processing_step:
                with open(channel) as r:
                    channel_dict = json.load(r)
                (
                    experiment,
                    period,
                    run,
                    datatype,
                    timestamp,
                    channel_name,
                    name,
                ) = os.path.basename(channel).split("-")
                out_dict[channel_name] = channel_dict
            else:
                pass

        pathlib.Path(os.path.dirname(out_file)).mkdir(parents=True, exist_ok=True)
        with open(out_file, "w") as w:
            json.dump(out_dict, w, indent=4)

    elif file_extension == ".pkl":
        out_dict = {}
        for channel in channel_files:
            if os.path.splitext(channel)[0].split("-")[-1] == processing_step:
                with open(channel, "rb") as r:
                    channel_dict = pkl.load(r)
                (
                    experiment,
                    period,
                    run,
                    datatype,
                    timestamp,
                    channel_name,
                    name,
                ) = os.path.basename(channel).split("-")
                out_dict[channel_name] = channel_dict
            else:
                pass
        pathlib.Path(os.path.dirname(out_file)).mkdir(parents=True, exist_ok=True)
        with open(out_file, "wb") as w:
            pkl.dump(out_dict, w, protocol=pkl.HIGHEST_PROTOCOL)

    elif file_extension == ".dat" or file_extension == ".dir":
        _out_file = os.path.splitext(out_file)[0]
        pathlib.Path(os.path.dirname(_out_file)).mkdir(parents=True, exist_ok=True)
        common_dict = {}
        with shelve.open(_out_file, "c", protocol=pkl.HIGHEST_PROTOCOL) as shelf:
            for channel in channel_files:
                if os.path.splitext(channel)[0].split("-")[-1] == processing_step:
                    with open(channel, "rb") as r:
                        channel_dict = pkl.load(r)
                    (
                        experiment,
                        period,
                        run,
                        datatype,
                        timestamp,
                        channel_name,
                        name,
                    ) = os.path.basename(channel).split("-")
                    if "common" in list(channel_dict):
                        chan_common_dict = channel_dict.pop("common")
                        common_dict[channel_name] = chan_common_dict
                    shelf[channel_name] = channel_dict
                else:
                    pass
            if len(common_dict) > 0:
                shelf["common"] = common_dict
