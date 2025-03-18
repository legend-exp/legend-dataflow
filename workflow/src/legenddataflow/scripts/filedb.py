import argparse
import logging
from pathlib import Path

import numpy as np
from dbetto.catalog import Props
from lgdo import lh5
from pygama.flow.file_db import FileDB


def build_filedb() -> None:
    argparser = argparse.ArgumentParser()
    argparser.add_argument("--config", required=True)
    argparser.add_argument("--scan-path", required=True)
    argparser.add_argument("--output", required=True)
    argparser.add_argument("--log")
    argparser.add_argument("--assume-nonsparse", action="store_true")
    args = argparser.parse_args()

    config = Props.read_from(args.config)

    if args.log is not None:
        Path(args.log).parent.mkdir(parents=True, exist_ok=True)
        logging.basicConfig(level=logging.DEBUG, filename=args.log, filemode="w")
    else:
        logging.basicConfig(level=logging.DEBUG)

    logging.getLogger("legendmeta").setLevel(logging.INFO)
    logging.getLogger("numba").setLevel(logging.INFO)
    logging.getLogger("parse").setLevel(logging.INFO)
    logging.getLogger("lgdo").setLevel(logging.INFO)
    logging.getLogger("h5py._conv").setLevel(logging.INFO)

    log = logging.getLogger(__name__)

    fdb = FileDB(config, scan=False)
    fdb.scan_files([args.scan_path])
    fdb.scan_tables_columns(dir_files_conform=True)

    # augment dataframe with earliest timestamp found in file

    default = np.finfo("float64").max
    timestamps = np.zeros(len(fdb.df), dtype="float64")

    for i, row in enumerate(fdb.df.itertuples()):
        store = lh5.LH5Store(
            base_path=f"{fdb.data_dir}/{fdb.tier_dirs['raw']}", keep_open=True
        )

        # list of first timestamps for each channel
        loc_timestamps = np.full(
            len(row.raw_tables), fill_value=default, dtype="float64"
        )

        msg = f"finding first timestamp in {fdb.data_dir}/{fdb.tier_dirs['raw']}/{row.raw_file}"
        log.info(msg)

        found = False
        for j, table in enumerate(row.raw_tables):
            try:
                loc_timestamps[j] = store.read(
                    fdb.table_format["raw"].format(ch=table) + "/timestamp",
                    row.raw_file.strip("/"),
                    n_rows=1,
                )[0]
                found = True
            except KeyError:
                pass

            if found and args.assume_nonsparse:
                break

        if (loc_timestamps == default).all() or not found:
            msg = "something went wrong! no valid first timestamp found. Likely: the file is empty"
            raise RuntimeError(msg)

        timestamps[i] = np.min(loc_timestamps)

        msg = f"found {timestamps[i]}"
        log.info(msg)

        if timestamps[i] < 0 or timestamps[i] > 4102444800:
            msg = f"something went wrong! timestamp {timestamps[i]} does not make sense"
            raise RuntimeError(msg)

    fdb.df["first_timestamp"] = timestamps

    fdb.to_disk(args.output, wo_mode="of")
