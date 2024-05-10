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
# TODO: optimize buffer_len
for chunk, entry, n_rows in lh5.LH5Iterator(args.evt_file, "evt", buffer_len=32768):
    evt = chunk.view_as("ak")

    # remove unwanted events
    skm = eval(f"evt[{evt_filter}]")
    # make it rectangular and make an LGDO Table
    out_table = Table(ak.flatten(skm, axis=1))

    # remove unwanted columns
    out_table_skm = Table(size=len(out_table))
    for field in out_fields:
        # table nesting is labeled by '.' in the config
        items = field.split(".")

        # add column to output table
        if field in out_fields:
            # get to actual nested field recursively
            ptr1 = out_table
            ptr2 = out_table_skm
            for item in items[:-1]:
                # make intermediate tables in new table
                if item not in ptr2:
                    ptr2.add_field(item, Table(size=(out_table)))

                # get non-table LGDO recursively
                ptr1 = ptr1[item]
                ptr2 = ptr2[item]

            # finally add column to new table
            ptr2.add_field(items[-1], ptr1[item])

            # forward LGDO attributes
            attrs = evt.flatten()[field.replace(".", "_")].attrs
            for attr, val in attrs:
                if attr != "datatype":
                    ptr2.attrs[attr] = val

    # write-append to disk
    store.write(
        out_table_skm, "skm", args.output, n_rows=n_rows, wo_mode=("w" if entry == 0 else "a")
    )
