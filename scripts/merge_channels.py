import os,json
import snakemake as smk
#from util.utils import *
import pathlib
import pickle as pkl

channel_files = snakemake.input

if isinstance(channel_files,str):
    if channel_files.split('.')[-1] == 'chanlist':
        with open(channel_files) as f:  
            channel_files = f.read().splitlines()
    else:
        channel_files = [channel_files]

#Combine par files
out_dict = {}
for channel in channel_files:
    with open(channel,"r") as r:
        channel_dict = json.load(r)
    experiment, period, run,datatype, timestamp,channel_name, name = os.path.basename(channel).split("-")

    out_dict[channel_name] = channel_dict

pathlib.Path(os.path.dirname(snakemake.output[0])).mkdir(parents=True, exist_ok=True)
with open(snakemake.output[0],"w") as w:
    json.dump(out_dict, w,indent=4)

###combine plot files

#convert pars files to plot files



#merge

#for channel in channel_files:
#    with open(channel,"rb") as r:
#        channel_dict = pkl.load(r)
#    experiment, period, run,datatype, timestamp,channel_name, name = os.path.basename(channel).split("-")

#    out_dict[channel_name] = channel_dict

#pathlib.Path(os.path.dirname(snakemake.output[1])).mkdir(parents=True, exist_ok=True)
#with open(snakemake.output[1],"wb") as w:
#    pkl.dump(plot_dict, w)
