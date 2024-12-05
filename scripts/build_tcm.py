import argparse
import logging
import logging.config
from pathlib import Path

import lgdo.lh5 as lh5
import numpy as np
from daq2lh5.orca import orca_flashcam
from legendmeta import TextDB
from legendmeta.catalog import Props
from pygama.evt.build_tcm import build_tcm

argparser = argparse.ArgumentParser()
argparser.add_argument("input", help="input file", type=str)
argparser.add_argument("output", help="output file", type=str)
argparser.add_argument("--datatype", help="Datatype", type=str, required=True)
argparser.add_argument("--timestamp", help="Timestamp", type=str, required=True)
argparser.add_argument("--configs", help="config file", type=str)
argparser.add_argument("--log", help="log file", type=str)
args = argparser.parse_args()

configs = TextDB(args.configs, lazy=True).on(args.timestamp, system=args.datatype)
config_dict = configs["snakemake_rules"]["tier_tcm"]
log_config = config_dict["options"]["logging"]

Path(args.log).parent.mkdir(parents=True, exist_ok=True)
log_config = Props.read_from(log_config)
log_config["handlers"]["file"]["filename"] = args.log
logging.config.dictConfig(log_config)
log = logging.getLogger("test")

Path(args.output).parent.mkdir(parents=True, exist_ok=True)


settings = Props.read_from(config_dict["inputs"]["config"])

rng = np.random.default_rng()
temp_output = f"{args.output}.{rng.integers(0, 99999):05d}"

# get the list of channels by fcid
ch_list = lh5.ls(args.input, "/ch*")
fcid_channels = {}
for ch in ch_list:
    key = int(ch[2:])
    fcid = orca_flashcam.get_fcid(key)
    if fcid not in fcid_channels:
        fcid_channels[fcid] = []
    fcid_channels[fcid].append(f"/{ch}/raw")

# make a hardware_tcm_[fcid] for each fcid
for fcid in fcid_channels:
    build_tcm(
        [(args.input, fcid_channels[fcid])],
        out_file=temp_output,
        out_name=f"hardware_tcm_{fcid}",
        wo_mode="o",
        **settings,
    )

Path(temp_output).rename(args.output)
