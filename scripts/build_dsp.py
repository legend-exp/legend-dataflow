import argparse
import re
import time
import warnings
from pathlib import Path

import numpy as np
from dspeed import build_dsp
from legendmeta import LegendMetadata, TextDB
from legendmeta.catalog import Props
from lgdo import lh5
from util.log import build_log


def replace_list_with_array(dic):
    for key, value in dic.items():
        if isinstance(value, dict):
            dic[key] = replace_list_with_array(value)
        elif isinstance(value, list):
            dic[key] = np.array(value, dtype="float32")
        else:
            pass
    return dic


warnings.filterwarnings(action="ignore", category=RuntimeWarning)

argparser = argparse.ArgumentParser()
argparser.add_argument("--configs", help="configs path", type=str, required=True)
argparser.add_argument("--metadata", help="metadata", type=str, required=True)
argparser.add_argument("--log", help="log file", type=str)

argparser.add_argument("--datatype", help="Datatype", type=str, required=True)
argparser.add_argument("--timestamp", help="Timestamp", type=str, required=True)
argparser.add_argument("--tier", help="Tier", type=str, required=True)

argparser.add_argument("--pars_file", help="database file for detector", nargs="*", default=[])
argparser.add_argument("--input", help="input file", type=str)

argparser.add_argument("--output", help="output file", type=str)
argparser.add_argument("--db_file", help="db file", type=str)
args = argparser.parse_args()

configs = TextDB(args.configs, lazy=True)
config_dict = configs.on(args.timestamp, system=args.datatype)["snakemake_rules"]
if args.tier in ["dsp", "psp"]:
    config_dict = config_dict["tier_dsp"]
elif args.tier in ["ann", "pan"]:
    config_dict = config_dict["tier_ann"]
else:
    msg = f"Tier {args.tier} not supported"
    raise ValueError(msg)

log = build_log(config_dict, args.log)

channel_dict = config_dict["inputs"]["processing_chain"]
settings_dict = config_dict["options"].get("settings", {})
if isinstance(settings_dict, str):
    settings_dict = Props.read_from(settings_dict)

meta = LegendMetadata(path=args.metadata)
chan_map = meta.channelmap(args.timestamp, system=args.datatype)

if args.tier in ["ann", "pan"]:
    channel_dict = {
        f"ch{chan_map[chan].daq.rawid:07}/dsp": Props.read_from(file)
        for chan, file in channel_dict.items()
    }
else:
    channel_dict = {
        f"ch{chan_map[chan].daq.rawid:07}/raw": Props.read_from(file)
        for chan, file in channel_dict.items()
    }
db_files = [
    par_file for par_file in args.pars_file if Path(par_file).suffix in (".json", ".yaml", ".yml")
]

database_dic = Props.read_from(db_files, subst_pathvar=True)
database_dic = replace_list_with_array(database_dic)

Path(args.output).parent.mkdir(parents=True, exist_ok=True)

rng = np.random.default_rng()
rand_num = f"{rng.integers(0, 99999):05d}"
temp_output = f"{args.output}.{rand_num}"

start = time.time()

build_dsp(
    args.input,
    temp_output,
    {},
    database=database_dic,
    chan_config=channel_dict,
    write_mode="r",
    buffer_len=settings_dict.get("buffer_len", 1000),
    block_width=settings_dict.get("block_width", 16),
)

log.info(f"build_dsp finished in {time.time()-start}")
Path(temp_output).rename(args.output)

key = Path(args.output).name.replace(f"-tier_{args.tier}.lh5", "")

if args.tier in ["dsp", "psp"]:

    raw_channels = [channel for channel in lh5.ls(args.input) if re.match("(ch\\d{7})", channel)]
    raw_fields = [field.split("/")[-1] for field in lh5.ls(args.input, f"{raw_channels[0]}/raw/")]

    outputs = {}
    channels = []
    for channel, chan_dict in channel_dict.items():
        output = chan_dict["outputs"]
        in_dict = False
        for entry in outputs:
            if outputs[entry]["fields"] == output:
                outputs[entry]["channels"].append(channel.split("/")[0])
                in_dict = True
        if in_dict is False:
            outputs[f"group{len(list(outputs))+1}"] = {
                "channels": [channel.split("/")[0]],
                "fields": output,
            }
        channels.append(channel.split("/")[0])

    full_dict = {
        "valid_fields": {
            "raw": {"group1": {"fields": raw_fields, "channels": raw_channels}},
            "dsp": outputs,
        },
        "valid_keys": {key: {"valid_channels": {"raw": raw_channels, "dsp": channels}}},
    }
else:
    outputs = {}
    channels = []
    for channel, chan_dict in channel_dict.items():
        output = chan_dict["outputs"]
        in_dict = False
        for entry in outputs:
            if outputs[entry]["fields"] == output:
                outputs[entry]["channels"].append(channel.split("/")[0])
                in_dict = True
        if in_dict is False:
            outputs[f"group{len(list(outputs))+1}"] = {
                "channels": [channel.split("/")[0]],
                "fields": output,
            }
        channels.append(channel.split("/")[0])

    full_dict = {
        "valid_fields": {
            "ann": outputs,
        },
        "valid_keys": {key: {"valid_channels": {"ann": channels}}},
    }

Path(args.db_file).parent.mkdir(parents=True, exist_ok=True)
Props.write_to(args.db_file, full_dict)
