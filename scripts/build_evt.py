import argparse
import json
import logging
import os
import pathlib
import time

import lgdo.lh5_store as lh5
from legendmeta import LegendMetadata
from legendmeta.catalog import Props
from pygama.evt.build_evt import build_evt

argparser = argparse.ArgumentParser()
argparser.add_argument("--hit_file", help="hit file", type=str)
argparser.add_argument("--dsp_file", help="dsp file", type=str)
argparser.add_argument("--tcm_file", help="tcm file", type=str)

argparser.add_argument("--configs", help="configs", type=str, required=True)
argparser.add_argument("--datatype", help="Datatype", type=str, required=True)
argparser.add_argument("--timestamp", help="Timestamp", type=str, required=True)
argparser.add_argument("--tier", help="Tier", type=str, required=True)

argparser.add_argument("--metadata", help="metadata path", type=str, required=True)

argparser.add_argument("--log", help="log_file", type=str)

argparser.add_argument("--output", help="output file", type=str)
args = argparser.parse_args()

pathlib.Path(os.path.dirname(args.log)).mkdir(parents=True, exist_ok=True)
logging.basicConfig(level=logging.DEBUG, filename=args.log, filemode="w")
logging.getLogger("numba").setLevel(logging.INFO)
logging.getLogger("parse").setLevel(logging.INFO)
logging.getLogger("lgdo").setLevel(logging.INFO)
logging.getLogger("h5py._conv").setLevel(logging.INFO)


log = logging.getLogger(__name__)

# load in config
configs = LegendMetadata(path=args.configs)
if args.tier == "evt" or args.tier == "pet":
    evt_config_file = configs.on(args.timestamp,  system=args.datatype)["snakemake_rules"]["tier_evt"][
        "inputs"
    ]["evt_config"]
else:
    msg = "unknown tier"
    raise ValueError(msg)

evt_config = Props.read_from(evt_config_file)

meta = LegendMetadata(path=args.metadata)
chmap = meta.channelmap(args.timestamp)

# block for snakemake to fill in channel lists
for field, dic in evt_config["channels"].items():
    if isinstance(dic, dict):
        chans = chmap.map("system", unique=False)[dic["system"]]
        if "selectors" in dic.keys():
            for k,val in dic["selectors"].items():
                chans = chans.map(k, unique=False)[val]
        chans = [f"ch{chan}" for chan in list(chans.map("daq.rawid"))]
        evt_config["channels"][field] = chans

log.debug(json.dumps(evt_config["channels"],indent=2))
log.debug(json.dumps(evt_config,indent=2))

t_start = time.time()
pathlib.Path(os.path.dirname(args.output)).mkdir(parents=True, exist_ok=True)
build_evt(f_tcm = args.tcm_file, 
        f_dsp = args.dsp_file, 
        f_hit = args.hit_file, 
        f_evt = args.output, 
        evt_config = evt_config,
        evt_group = "evt",
        tcm_group = "hardware_tcm_1",
        dsp_group = "dsp",
        hit_group = "hit",
        tcm_id_table_pattern = "ch{}",
        wo_mode="o")
t_elap = time.time() - t_start
log.info(f"Done!  Time elapsed: {t_elap:.2f} sec.")

