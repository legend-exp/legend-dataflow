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

if args.log is not None:
    pathlib.Path(os.path.dirname(args.log)).mkdir(parents=True, exist_ok=True)
    logging.basicConfig(level=logging.DEBUG, filename=args.log, filemode="w")
else:
    logging.basicConfig(level=logging.DEBUG)

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

            evt_config[key] = _evt_config
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
    datainfo = {
        "tcm": (args.tcm_file, "hardware_tcm_1", "ch{}"),
        "dsp": (args.dsp_file, "dsp", "ch{}"),
        "hit": (args.hit_file, "hit", "ch{}"),
        "evt": (None, "evt"),
    }

    tables[key] = build_evt(
        datainfo,
        config,
    )

tbl = Table(col_dict=tables)
sto.write(obj=tbl, name="evt", lh5_file=temp_output, wo_mode="a")

os.rename(temp_output, args.output)
t_elap = time.time() - t_start
log.info(f"Done!  Time elapsed: {t_elap:.2f} sec.")
