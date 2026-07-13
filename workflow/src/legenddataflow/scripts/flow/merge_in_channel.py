"""Console script ``merge-in-channels``: merge per-channel parameter outputs
into a single channel-keyed output file."""

from __future__ import annotations

import argparse
import pickle as pkl
from pathlib import Path

import lh5
from dbetto import AttrsDict
from dbetto.catalog import Props
from legenddataflowscripts.utils import check_input_files
from lgdo.types import Struct

from legenddataflow.methods import ChannelProcKey

from .utils import replace_path


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

    if len(input_files) == 0:
        msg = "no input files provided, cannot merge anything"
        raise ValueError(msg)
    check_input_files(input_files, "--input")

    if args.out_db and not args.in_db:
        msg = "--out-db requires --in-db: there is no input database to update"
        raise ValueError(msg)

    file_extension = Path(args.output).suffix
    out_file = args.output

    Path(args.output).parent.mkdir(parents=True, exist_ok=True)

    if file_extension in (".json", ".yaml", ".yml"):
        Props.write_to(out_file, AttrsDict(Props.read_from(input_files)).to_dict())

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
        channel = None
        for infile in input_files:
            fkey = ChannelProcKey.get_filekey_from_pattern(Path(infile).name)
            if fkey.identifier is None:
                msg = (
                    f"cannot derive the output group for {infile}: processing "
                    f"step {fkey.processing_step!r} has no identifier part "
                    "(expected {type}_{tier}_{identifier})"
                )
                raise ValueError(msg)
            if channel is not None and fkey.channel != channel:
                msg = (
                    f"input files span multiple channels ({channel} and "
                    f"{fkey.channel}): merge-in-channels merges processing "
                    "steps of a single channel"
                )
                raise ValueError(msg)
            channel = fkey.channel

            tb_in = lh5.read(fkey.channel, infile)
            if fkey.identifier not in tb_in:
                msg = (
                    f"table {fkey.channel!r} in {infile} has no field "
                    f"{fkey.identifier!r}: available fields are {list(tb_in)}"
                )
                raise KeyError(msg)
            objects[fkey.identifier] = tb_in[fkey.identifier]

        lh5.write(
            Struct(objects),
            name=channel,
            lh5_file=out_file,
            wo_mode="a",
        )

        if args.out_db:
            for infile in input_files:
                db_dict = replace_path(db_dict, infile, args.output)
            Props.write_to(args.out_db, AttrsDict(db_dict).to_dict())

    else:
        msg = (
            f"unsupported output file extension {file_extension!r} for {args.output}: "
            "expected one of .json, .yaml, .yml, .pkl or .lh5"
        )
        raise ValueError(msg)
