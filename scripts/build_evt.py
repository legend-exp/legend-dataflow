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
from lgdo.types import Array
from pygama.evt import build_evt

sto = lh5.LH5Store()


def find_matching_values_with_delay(arr1, arr2, jit_delay):
    matching_values = []

    # Create an array with all possible delay values
    delays = np.arange(0, int(1e9 * jit_delay)) * jit_delay

    for delay in delays:
        arr2_delayed = arr2 + delay

        # Find matching values and indices
        mask = np.isin(arr1, arr2_delayed, assume_unique=True)
        matching_values.extend(arr1[mask])

    return np.unique(matching_values)


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

logging.getLogger("legendmeta").setLevel(logging.INFO)
logging.getLogger("numba").setLevel(logging.INFO)
logging.getLogger("parse").setLevel(logging.INFO)
logging.getLogger("lgdo").setLevel(logging.INFO)
logging.getLogger("h5py._conv").setLevel(logging.INFO)


log = logging.getLogger(__name__)

# load in config
configs = LegendMetadata(path=args.configs)
if args.tier == "evt" or args.tier == "pet":
    config_dict = configs.on(args.timestamp, system=args.datatype)["snakemake_rules"]["tier_evt"][
        "inputs"
    ]
    evt_config_file = config_dict["evt_config"]
else:
    msg = "unknown tier"
    raise ValueError(msg)

meta = LegendMetadata(path=args.metadata)
chmap = meta.channelmap(args.timestamp)

evt_config = Props.read_from(evt_config_file)

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

log.debug(json.dumps(evt_config["channels"], indent=2))

t_start = time.time()
pathlib.Path(os.path.dirname(args.output)).mkdir(parents=True, exist_ok=True)

rng = np.random.default_rng()
rand_num = f"{rng.integers(0,99999):05d}"
temp_output = f"{args.output}.{rand_num}"

table = build_evt(
    {
        "tcm": (args.tcm_file, "hardware_tcm_1", "ch{}"),
        "dsp": (args.dsp_file, "dsp", "ch{}"),
        "hit": (args.hit_file, "hit", "ch{}"),
        "evt": (temp_output, "evt"),
    },
    evt_config,
)

if "muon_config" in config_dict and config_dict["muon_config"] is not None:
    muon_config = Props.read_from(config_dict["muon_config"]["evt_config"])
    field_config = Props.read_from(config_dict["muon_config"]["field_config"])
    # block for snakemake to fill in channel lists
    for field, dic in muon_config["channels"].items():
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
            muon_config["channels"][field] = chans

    trigger_timestamp = table[field_config["ged_timestamp"]["table"]][
        field_config["ged_timestamp"]["field"]
    ].nda
    if "hardware_tcm_2" in lh5.ls(args.tcm_file):
        muon_table = build_evt(
            {
                "tcm": (args.tcm_file, "hardware_tcm_2", "ch{}"),
                "dsp": (args.dsp_file, "dsp", "ch{}"),
                "hit": (args.hit_file, "hit", "ch{}"),
                "evt": (None, "evt"),
            },
            muon_config,
        )

        muon_timestamp = muon_table[field_config["muon_timestamp"]["field"]].nda
        muon_tbl_flag = muon_table[field_config["muon_flag"]["field"]].nda
        if len(muon_timestamp[muon_tbl_flag]) > 0:
            is_muon_veto_triggered = find_matching_values_with_delay(
                trigger_timestamp, muon_timestamp[muon_tbl_flag], field_config["jitter"]
            )
            muon_flag = np.isin(trigger_timestamp, is_muon_veto_triggered)
        else:
            muon_flag = np.zeros(len(trigger_timestamp), dtype=bool)
    else:
        muon_flag = np.zeros(len(trigger_timestamp), dtype=bool)
    table[field_config["output_field"]["table"]].add_column(
        field_config["output_field"]["field"], Array(muon_flag)
    )

sto.write(obj=table, name="evt", lh5_file=temp_output, wo_mode="a")

os.rename(temp_output, args.output)
t_elap = time.time() - t_start
log.info(f"Done!  Time elapsed: {t_elap:.2f} sec.")
