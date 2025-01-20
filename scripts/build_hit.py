import argparse
import time
from pathlib import Path

from legendmeta import LegendMetadata, TextDB
from legendmeta.catalog import Props
from lgdo import lh5
from pygama.hit.build_hit import build_hit
from util.log import build_log

argparser = argparse.ArgumentParser()
argparser.add_argument("--input", help="input file", type=str)
argparser.add_argument("--pars_file", help="hit pars file", nargs="*")

argparser.add_argument("--configs", help="configs", type=str, required=True)
argparser.add_argument("--metadata", help="metadata", type=str, required=True)
argparser.add_argument("--log", help="log_file", type=str)

argparser.add_argument("--datatype", help="Datatype", type=str, required=True)
argparser.add_argument("--timestamp", help="Timestamp", type=str, required=True)
argparser.add_argument("--tier", help="Tier", type=str, required=True)

argparser.add_argument("--output", help="output file", type=str)
argparser.add_argument("--db_file", help="db file", type=str)
args = argparser.parse_args()

configs = TextDB(args.configs, lazy=True)
if args.tier == "hit" or args.tier == "pht":
    config_dict = configs.on(args.timestamp, system=args.datatype)["snakemake_rules"]["tier_hit"]
else:
    msg = "unknown tier"
    raise ValueError(msg)

log = build_log(config_dict, args.log)

channel_dict = config_dict["inputs"]["hit_config"]
settings_dict = config_dict["options"].get("settings", {})
if isinstance(settings_dict, str):
    settings_dict = Props.read_from(settings_dict)

meta = LegendMetadata(path=args.metadata)
chan_map = meta.channelmap(args.timestamp, system=args.datatype)

pars_dict = Props.read_from(args.pars_file)
pars_dict = {chan: chan_dict["pars"] for chan, chan_dict in pars_dict.items()}

hit_dict = {}
channels_present = lh5.ls(args.input)
for channel in pars_dict:
    chan_pars = pars_dict[channel].copy()
    try:
        detector = chan_map.map("daq.rawid")[int(channel[2:])].name
        if detector in channel_dict:
            cfg_dict = Props.read_from(channel_dict[detector])
            Props.add_to(cfg_dict, chan_pars)
            chan_pars = cfg_dict

        if channel in channels_present:
            hit_dict[f"{channel}/dsp"] = chan_pars
    except KeyError:
        pass

t_start = time.time()
Path(args.output).parent.mkdir(parents=True, exist_ok=True)
build_hit(args.input, lh5_tables_config=hit_dict, outfile=args.output)
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

key = args.output.replace(f"-tier_{args.tier}.lh5", "")

full_dict = {
    "valid_fields": {args.tier: hit_outputs},
    "valid_keys": {key: {"valid_channels": {args.tier: hit_channels}}},
}

Path(args.db_file).parent.mkdir(parents=True, exist_ok=True)
Props.write_to(args.db_file, full_dict)
