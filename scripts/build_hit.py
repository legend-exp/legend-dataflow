import argparse, pathlib
import numpy as np
import os,json

from legendmeta import LegendMetadata
from legendmeta.catalog import Props

from pygama.hit.build_hit import build_hit
import lgdo.lh5_store as lh5

import time
import logging

argparser = argparse.ArgumentParser()
argparser.add_argument("--input", help="input file", type=str)
argparser.add_argument("--pars_file", help="hit pars file", nargs="*")

argparser.add_argument("--configs", help="configs", type=str, required=True)
argparser.add_argument("--datatype", help="Datatype", type=str, required=True)
argparser.add_argument("--timestamp", help="Timestamp", type=str, required=True)
argparser.add_argument("--tier", help="Tier", type=str, required=True)

argparser.add_argument("--log", help="log_file", type=str)

argparser.add_argument("--output", help="output file", type=str)
argparser.add_argument("--db_file", help="db file", type=str)
args = argparser.parse_args()

pathlib.Path(os.path.dirname(args.log)).mkdir(parents=True, exist_ok=True)
logging.basicConfig(level=logging.DEBUG, filename=args.log, filemode='w')
logging.getLogger('numba').setLevel(logging.INFO)
logging.getLogger('parse').setLevel(logging.INFO)
logging.getLogger('pygama.lgdo.lh5_store').setLevel(logging.INFO)
logging.getLogger('h5py._conv').setLevel(logging.INFO)


log = logging.getLogger(__name__)


configs = LegendMetadata(path = args.configs)
if args.tier == "hit":  
    channel_dict = configs.on(args.timestamp, system=args.datatype)['snakemake_rules']['tier_hit']["inputs"]['hit_config']
if args.tier == "pht": 
    channel_dict = configs.on(args.timestamp, system=args.datatype)['snakemake_rules']['tier_hit']["inputs"]['hit_config']
else:
    raise ValueError("unknown tier")

if isinstance(args.pars_file, list):
    pars_dict = Props.read_from(args.pars_file)
else:
    with open(args.pars_file) as f:
        pars_dict = json.load(f)

pars_dict = {chan: chan_dict["pars"] for chan, chan_dict in pars_dict.items()}

hit_dict ={}
channels_present = lh5.ls(args.input)
for channel in pars_dict:
    if channel in channel_dict:
        with open(channel_dict[channel], "r") as r:
            cfg_dict = json.load(r)
        Props.add_to(pars_dict[channel], cfg_dict)
    if channel in channels_present:
        hit_dict[f'{channel}/dsp']=pars_dict[channel]

t_start = time.time()
pathlib.Path(os.path.dirname(args.output)).mkdir(parents=True, exist_ok=True)
build_hit(args.input, lh5_tables_config=hit_dict, outfile =args.output)
t_elap = time.time() - t_start
log.info(f'Done!  Time elapsed: {t_elap:.2f} sec.')

hit_outputs = {}
hit_channels=[]
for channel,file in channel_dict.items():
    output = Props.read_from(file)["outputs"]
    in_dict = False
    for entry in hit_outputs:
        if hit_outputs[entry]["fields"]==output:
            hit_outputs[entry]["channels"].append(channel)
            in_dict=True
    if in_dict == False:
        hit_outputs[f"group{len(list(hit_outputs))+1}"]={"channels":[channel],
                                            "fields":output}
    hit_channels.append(channel)

key = os.path.basename(args.output).replace(f"-tier_{args.tier}.lh5","")

full_dict = {"valid_fields":{
    args.tier:hit_outputs
},
             
    "valid_keys":{
        key:{
            "valid_channels":{
                args.tier:hit_channels
            }
        }
    }
}

pathlib.Path(os.path.dirname(args.db_file)).mkdir(parents=True, exist_ok=True)
with open(args.db_file ,"w") as w:
    json.dump(full_dict, w, indent=4)