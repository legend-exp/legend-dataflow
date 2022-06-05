import os,json
import snakemake as smk
from utils import *
import pathlib

channel_files = snakemake.input

if isinstance(channel_files,str):
    if channel_files.split('.')[-1] == 'chanlist':
        with open(channel_files) as f:  
            channel_files = f.read().splitlines()
    else:
        channel_files = [channel_files]

out_dict = {}
for channel in channel_files:
    with open(channel,"r") as r:
        channel_dict = json.load(r)
    name, experiment, period, run,channel_name, tier = os.path.basename(channel).split("-")

    out_dict[channel_name] = channel_dict

pathlib.Path(os.path.dirname(snakemake.output[0])).mkdir(parents=True, exist_ok=True)
with open(snakemake.output[0],"w") as w:
    json.dump(out_dict, w,indent=4)