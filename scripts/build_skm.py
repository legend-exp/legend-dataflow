import argparse
import logging
import os
import pathlib

import awkward as ak
from legendmeta import TextDB
from legendmeta.catalog import Props
from lgdo import lh5
from lgdo.types import Table

argparser = argparse.ArgumentParser()
argparser.add_argument("--evt-file", help="evt file", required=True)
argparser.add_argument("--configs", help="configs", required=True)
argparser.add_argument("--datatype", help="datatype", required=True)
argparser.add_argument("--timestamp", help="timestamp", required=True)
argparser.add_argument("--log", help="log file", default=None)
argparser.add_argument("--output", help="output file", required=True)
args = argparser.parse_args()

if args.log is not None:
    pathlib.Path(os.path.dirname(args.log)).mkdir(parents=True, exist_ok=True)

logging.basicConfig(level=logging.DEBUG, filename=args.log, filemode="w")

logging.getLogger("numba").setLevel(logging.INFO)
logging.getLogger("parse").setLevel(logging.INFO)
logging.getLogger("lgdo").setLevel(logging.INFO)
logging.getLogger("h5py._conv").setLevel(logging.INFO)

log = logging.getLogger(__name__)

# load in config
configs = TextDB(args.configs, lazy=True).on(args.timestamp, system=args.datatype)
skm_config_file = configs["snakemake_rules"]["tier_skm"]["inputs"]["skm_config"]

evt_filter = Props.read_from(skm_config_file)["evt_filter"]
out_fields = Props.read_from(skm_config_file)["keep_fields"]

store = lh5.LH5Store(keep_open=True)

# read concatenated evt file in chunks to keep memory usage low
for chunk, entry, n_rows in lh5.LH5Iterator(args.evt_file, "evt", buffer_len=32768):
    evt = chunk.view_as("ak")

    # remove unwanted events
    skm = eval(f"evt[{evt_filter}]")
    rect_skm = ak.flatten(skm, axis=-1)

    # remove unwanted columns
    table = Table(rect_skm)
    for field in table:
        if field not in out_fields:
            ptr = table
            items = field.split(".")
            for item in items[:-1]:
                ptr = ptr[item]
            ptr.remove_column(items[-1])

    # write-append to disk
    store.write(table, "skm", args.output, n_rows=n_rows, wo_mode=("w" if entry == 0 else "a"))
