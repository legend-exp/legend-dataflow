from __future__ import annotations

import argparse
import pickle as pkl
from pathlib import Path

from dbetto.catalog import Props
from lgdo import lh5
from lgdo.types import Struct

from legenddataflow.methods import ChannelProcKey
from .utils import replace_path


def replace_path(d, old_path, new_path):
    if isinstance(d, dict):
        for k, v in d.items():
            d[k] = replace_path(v, old_path, new_path)
    elif isinstance(d, list):
        for i in range(len(d)):
            d[i] = replace_path(d[i], old_path, new_path)
    elif isinstance(d, str) and old_path in d:
        d = d.replace(old_path, new_path)
        d = d.replace(new_path, f"$_/{Path(new_path).name}")
    return d


def merge_in_channel() -> None:
    argparser = argparse.ArgumentParser()
    argparser.add_argument(
        "--input", help="input file", nargs="*", type=str, required=True
    )
    argparser.add_argument("--output", help="output file", type=str, required=True)
    argparser.add_argument(
        "--in-db",
        help="in db file (used for when lh5 files referred to in db)",
        type=str,
        required=False,
    )
    argparser.add_argument(
        "--out-db",
        help="lh5 file (used for when lh5 files referred to in db)",
        type=str,
        required=False,
    )
    args = argparser.parse_args()

    # change to only have 1 output file for multiple inputs
    # don't care about processing step, check if extension matches

    input_files = args.input.infiles if hasattr(args.input, "infiles") else args.input

    file_extension = Path(args.output).suffix
    out_file = args.output

    Path(args.output).parent.mkdir(parents=True, exist_ok=True)

    if file_extension in (".json", ".yaml", ".yml"):
        Props.write_to(out_file, Props.read_from(input_files))

    elif file_extension == ".pkl":
        out_dict = {}
        for infile in input_files:
            with Path(infile).open("rb") as r:
                indict = pkl.load(r)
            out_dict.update(indict)

        with Path(out_file).open("wb") as w:
            pkl.dump(out_dict, w, protocol=pkl.HIGHEST_PROTOCOL)

    elif file_extension == ".lh5":
        if args.in_db:
            db_dict = Props.read_from(args.in_db)
        objects = {}
        for infile in input_files:
            fkey = ChannelProcKey.get_filekey_from_pattern(Path(infile).name)
            tb_in = lh5.read(f"{fkey.channel}", infile)

            objects[f"{fkey.processing_step}".split("_", maxsplit=3)[2]] = tb_in

        lh5.write(
            Struct(objects),
            name=fkey.channel,
            lh5_file=out_file,
            wo_mode="a",
        )

        if args.out_db:
            if args.in_db:
                db_dict = replace_path(
                    db_dict, infile, args.output
                )
            
            Props.write_to(args.out_db, db_dict)