import argparse
import json
import os
import pathlib
import pickle as pkl
import shelve

import lgdo.lh5_store as lh5
from lgdo import Array
sto = lh5.LH5Store()

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
                
                for key in channel_dict.keys():
                    key_dict = channel_dict[key]
                    for key_pars in key_dict.keys():
                        if isinstance(key_dict[key_pars], str) and ("loadlh5" in key_dict[key_pars]):
                            out_lh5 = out_file.replace(".json",".lh5")
                            out_dict[channel_name][key][key_pars] = f"loadlh5('{out_lh5}', '{channel_name}/{key}')"
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
                    if isinstance(channel_dict, dict) and "common" in list(channel_dict):
                        chan_common_dict = channel_dict.pop("common")
                        common_dict[channel_name] = chan_common_dict
                    shelf[channel_name] = channel_dict
                else:
                    pass
            if len(common_dict) > 0:
                shelf["common"] = common_dict

    elif file_extension == ".lh5":
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
                
                for key in channel_dict.keys():
                    key_dict = channel_dict[key]
                    for key_pars in key_dict.keys():
                        if isinstance(key_dict[key_pars], str) and ("loadlh5" in key_dict[key_pars]):
                            path_to_file = key_dict[key_pars].split("'")[1]
                            path_in_file = key_dict[key_pars].split("'")[3]
                            data = sto.read_object(path_in_file, path_to_file)[0].nda
                            sto.write_object(
                                Array(data),
                                name=key,
                                lh5_file=out_file,
                                wo_mode="overwrite",
                                group=channel_name
                            )
            else:
                pass