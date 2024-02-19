import argparse
import json
import logging
import os
import pathlib
import time

import lgdo.lh5 as lh5
import numpy as np
from legendmeta import LegendMetadata
from legendmeta.catalog import Props
from lgdo.types import Table
from pygama.evt.build_evt import build_evt

sto = lh5.LH5Store()


def replace_evt_with_key(dic, new_key):
    for key, d in dic.items():
        if isinstance(d, dict):
            dic[key] = replace_evt_with_key(d, new_key)
        elif isinstance(d, list):
            dic[key] = [item.replace("evt", new_key) for item in d]
        elif isinstance(d, str):
            dic[key] = d.replace("evt", new_key)
        else:
            pass
    return dic


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
    evt_config_file = configs.on(args.timestamp, system=args.datatype)["snakemake_rules"][
        "tier_evt"
    ]["inputs"]["evt_config"]
else:
    msg = "unknown tier"
    raise ValueError(msg)

meta = LegendMetadata(path=args.metadata)
chmap = meta.channelmap(args.timestamp)

if isinstance(evt_config_file, dict):
    evt_config = {}
    for key, _evt_config in evt_config_file.items():
        if _evt_config is not None:
            _evt_config = Props.read_from(_evt_config)
            # block for snakemake to fill in channel lists
            for field, dic in _evt_config["channels"].items():
                if isinstance(dic, dict):
                    chans = chmap.map("system", unique=False)[dic["system"]]
                    if "selectors" in dic:
                        try:
                            for k, val in dic["selectors"].items():
                                chans = chans.map(k, unique=False)[val]
                        except KeyError:
                            chans = None
                    if chans is not None:
                        chans = [f"ch{chan}" for chan in list(chans.map("daq.rawid"))]
                    else:
                        chans = []
                    _evt_config["channels"][field] = chans
            evt_config[key] = replace_evt_with_key(_evt_config, f"evt/{key}")
else:
    evt_config = {"all": Props.read_from(evt_config_file)}
    # block for snakemake to fill in channel lists
    for field, dic in evt_config["channels"].items():
        if isinstance(dic, dict):
            chans = chmap.map("system", unique=False)[dic["system"]]
            if "selectors" in dic:
                try:
                    for k, val in dic["selectors"].items():
                        chans = chans.map(k, unique=False)[val]
                except KeyError:
                    chans = None
            if chans is not None:
                chans = [f"ch{chan}" for chan in list(chans.map("daq.rawid"))]
            else:
                chans = []
            evt_config["channels"][field] = chans

log.debug(json.dumps(evt_config, indent=2))

t_start = time.time()
pathlib.Path(os.path.dirname(args.output)).mkdir(parents=True, exist_ok=True)

rng = np.random.default_rng()
rand_num = f"{rng.integers(0,99999):05d}"
temp_output = f"{args.output}.{rand_num}"

tables = {}
for key, config in evt_config.items():
    tables[key] = build_evt(
        f_tcm=args.tcm_file,
        f_dsp=args.dsp_file,
        f_hit=args.hit_file,
        f_evt=None,
        evt_config=config,
        evt_group=f"evt/{key}" if key != "all" else "evt",
        tcm_group="hardware_tcm_1",
        dsp_group="dsp",
        hit_group="hit",
        tcm_id_table_pattern="ch{}"
    )

tbl = Table(col_dict=tables)
sto.write(obj=tbl, name="evt", lh5_file=temp_output, wo_mode="a")

os.rename(temp_output, args.output)
t_elap = time.time() - t_start
log.info(f"Done!  Time elapsed: {t_elap:.2f} sec.")
