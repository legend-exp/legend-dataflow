import os,json
import pathlib
import pickle as pkl
import numpy as np
import argparse
import pygama

argparser = argparse.ArgumentParser()
argparser.add_argument("--input", help="input file", nargs='*',type=str)
argparser.add_argument("--output", help="output file", nargs='*',type=str)
args = argparser.parse_args()



channel_files = args.input
for i,out_file in enumerate(args.output):
    file_extension = pathlib.Path(out_file).suffix
    processing_step = out_file.split("-")[-1]
    if file_extension == ".json":
        out_dict = {}
        for channel in channel_files:
            if channel.split("-")[-1] == processing_step:
                with open(channel,"r") as r:
                    channel_dict = json.load(r)
                experiment, period, run,datatype, timestamp,channel_name, name = os.path.basename(channel).split("-")
                out_dict[channel_name] = channel_dict
            else:
                pass
            
        pathlib.Path(os.path.dirname(out_file)).mkdir(parents=True, exist_ok=True)
        with open(out_file,"w") as w:
            json.dump(out_dict, w,indent=4)

    elif file_extension == ".pkl":
        out_dict = {}
        for channel in channel_files:
            if channel.split("-")[-1] == processing_step:
                with open(channel,"rb") as r:
                    channel_dict = pkl.load(r)
                experiment, period, run,datatype, timestamp,channel_name, name = os.path.basename(channel).split("-")
                out_dict[channel_name] = channel_dict
            else:
                pass
        pathlib.Path(os.path.dirname(out_file)).mkdir(parents=True, exist_ok=True)
        with open(out_file,"wb") as w:
            pkl.dump(out_dict, w)