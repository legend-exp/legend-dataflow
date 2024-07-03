import argparse
import logging
import os
import pathlib
import time

from legendmeta import TextDB
from legendmeta.catalog import Props
from lgdo import lh5
from pygama.hit.build_hit import build_hit
from util.utils import as_ro

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
logging.basicConfig(level=logging.DEBUG, filename=args.log, filemode="w")
logging.getLogger("numba").setLevel(logging.INFO)
logging.getLogger("parse").setLevel(logging.INFO)
logging.getLogger("lgdo").setLevel(logging.INFO)
logging.getLogger("h5py._conv").setLevel(logging.INFO)

log = logging.getLogger(__name__)

configs = TextDB(as_ro(args.configs), lazy=True)
if args.tier == "hit" or args.tier == "pht":
    channel_dict = configs.on(args.timestamp, system=args.datatype)["snakemake_rules"]["tier_hit"][
        "inputs"
    ]["hit_config"]
else:
    msg = "unknown tier"
    raise ValueError(msg)

pars_dict = Props.read_from(as_ro(args.pars_file))

pars_dict = {chan: chan_dict["pars"] for chan, chan_dict in pars_dict.items()}

hit_dict = {}
channels_present = lh5.ls(as_ro(args.input))
for channel in pars_dict:
    chan_pars = pars_dict[channel].copy()
    if channel in channel_dict:
        cfg_dict = Props.read_from(channel_dict[channel])
        Props.add_to(cfg_dict, chan_pars)
        chan_pars = cfg_dict

    if channel in channels_present:
        hit_dict[f"{channel}/dsp"] = chan_pars

t_start = time.time()
pathlib.Path(os.path.dirname(args.output)).mkdir(parents=True, exist_ok=True)
build_hit(as_ro(args.input), lh5_tables_config=hit_dict, outfile=args.output)
t_elap = time.time() - t_start
log.info(f"Done!  Time elapsed: {t_elap:.2f} sec.")

hit_outputs = {}
hit_channels = []
for channel, file in channel_dict.items():
    output = Props.read_from(file)["outputs"]
    in_dict = False
    for entry in hit_outputs:
        if hit_outputs[entry]["fields"] == output:
            hit_outputs[entry]["channels"].append(channel)
            in_dict = True
    if in_dict is False:
        hit_outputs[f"group{len(list(hit_outputs))+1}"] = {
            "channels": [channel],
            "fields": output,
        }
    hit_channels.append(channel)

key = os.path.basename(args.output).replace(f"-tier_{args.tier}.lh5", "")

full_dict = {
    "valid_fields": {args.tier: hit_outputs},
    "valid_keys": {key: {"valid_channels": {args.tier: hit_channels}}},
}

pathlib.Path(os.path.dirname(args.db_file)).mkdir(parents=True, exist_ok=True)
Props.write_to(args.db_file, full_dict)
