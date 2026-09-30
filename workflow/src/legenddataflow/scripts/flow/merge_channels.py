"""Console script ``merge-channels``: merge per-channel parameter or plot
outputs for one timestamp into a single per-tier file (yaml/json,
pickle/shelve or lh5)."""

from __future__ import annotations

import argparse
import dbm.dumb
import pickle as pkl
import shelve
from pathlib import Path

import lh5
from dbetto import AttrsDict
from dbetto.catalog import Props
from legenddataflowscripts.utils import check_input_files

from legenddataflow.methods import ChannelProcKey

from .utils import replace_path


def merge_plot_data(channel_files, out_file) -> None:
    """Write the data (not the figures) of per-channel plot pickles to one lh5 file.

    Each channel lands under ``/<channel>``; ``common`` entries are gathered
    under ``/common/<channel>``, mirroring the plot shelf layout.
    """
    # only this mode needs plot_dict_to_lgdo (legend-dataflow-scripts >= 0.4.0)
    from legenddataflowscripts.utils import plot_dict_to_lgdo  # noqa: PLC0415
    from lgdo import Struct  # noqa: PLC0415

    common = {}
    for channel in channel_files:
        if Path(channel).suffix != ".pkl":
            msg = f"input file {channel} is not a plot pickle"
            raise RuntimeError(msg)
        with Path(channel).open("rb") as r:
            channel_dict = pkl.load(r)
        fkey = ChannelProcKey.get_filekey_from_pattern(Path(channel).name)
        if "common" in channel_dict:
            common[fkey.channel] = channel_dict.pop("common")
        struct = plot_dict_to_lgdo(channel_dict)
        if struct is not None:
            lh5.write(struct, name=fkey.channel, lh5_file=out_file, wo_mode="append")
    common_struct = plot_dict_to_lgdo(common)
    if common_struct is not None:
        lh5.write(common_struct, name="common", lh5_file=out_file, wo_mode="append")
    if not Path(out_file).exists():  # nothing convertible, still leave an output
        lh5.write(Struct({}), name="empty", lh5_file=out_file, wo_mode="append")


def merge_channels() -> None:
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

    channel_files = args.input.infiles if hasattr(args.input, "infiles") else args.input

    if len(channel_files) == 0:
        msg = "no input files provided, cannot merge anything"
        raise ValueError(msg)
    check_input_files(channel_files, "--input")

    if args.out_db and not args.in_db:
        msg = "--out-db requires --in-db: there is no input database to update"
        raise ValueError(msg)

    file_extension = Path(args.output).suffix

    if file_extension in (".dat", ".dir"):
        common_dict = {}
        out_file = Path(args.output).with_suffix("")
    else:
        out_file = args.output

    Path(args.output).parent.mkdir(parents=True, exist_ok=True)

    if file_extension in (".json", ".yaml", ".yml"):
        out_dict = {}
        for channel in channel_files:
            if Path(channel).suffix == file_extension:
                channel_dict = Props.read_from(channel)
                fkey = ChannelProcKey.get_filekey_from_pattern(Path(channel).name)
                out_dict[fkey.channel] = AttrsDict(channel_dict).to_dict()
            else:
                msg = (
                    f"input file {channel} does not match the extension "
                    f"{file_extension!r} of the output file {args.output}"
                )
                raise RuntimeError(msg)

        Props.write_to(out_file, out_dict)

    elif file_extension == ".pkl":
        out_dict = {}
        for channel in channel_files:
            with Path(channel).open("rb") as r:
                channel_dict = pkl.load(r)
            fkey = ChannelProcKey.get_filekey_from_pattern(Path(channel).name)
            out_dict[fkey.channel] = channel_dict

        with Path(out_file).open("wb") as w:
            pkl.dump(out_dict, w, protocol=pkl.HIGHEST_PROTOCOL)

    elif file_extension in (".dat", ".dir"):
        common_dict = {}
        # Open the 'dumb' database directly.
        # This forces a two-file backend, avoiding the automatic selection
        # of a one-file backend (e.g. xxx.db) which breaks snakemake's output
        # file detection
        with (
            dbm.dumb.open(str(out_file), "c") as db_object,
            shelve.Shelf(db_object, protocol=pkl.HIGHEST_PROTOCOL) as shelf,
        ):
            for channel in channel_files:
                with Path(channel).open("rb") as r:
                    channel_dict = pkl.load(r)
                fkey = ChannelProcKey.get_filekey_from_pattern(Path(channel).name)
                if isinstance(channel_dict, dict) and "common" in list(channel_dict):
                    chan_common_dict = channel_dict.pop("common")
                    common_dict[fkey.channel] = chan_common_dict
                shelf[fkey.channel] = channel_dict
            if len(common_dict) > 0:
                shelf["common"] = common_dict

    elif file_extension == ".lh5" and Path(channel_files[0]).suffix == ".pkl":
        merge_plot_data(channel_files, out_file)

    elif file_extension == ".lh5":
        if args.in_db:
            db_dict = Props.read_from(args.in_db)
        for channel in channel_files:
            if Path(channel).suffix == file_extension:
                fkey = ChannelProcKey.get_filekey_from_pattern(Path(channel).name)
                tb_in = lh5.read(f"{fkey.channel}", channel)

                lh5.write(
                    tb_in,
                    name=fkey.channel,
                    lh5_file=out_file,
                    wo_mode="a",
                )
                if args.in_db:
                    if fkey.channel not in db_dict:
                        msg = (
                            f"channel {fkey.channel} (from {channel}) not found in "
                            f"the input database {args.in_db}"
                        )
                        raise KeyError(msg)
                    db_dict[fkey.channel] = replace_path(
                        db_dict[fkey.channel], channel, args.output
                    )
            else:
                msg = (
                    f"input file {channel} does not match the extension "
                    f"{file_extension!r} of the output file {args.output}"
                )
                raise RuntimeError(msg)
        if args.out_db:
            Props.write_to(args.out_db, AttrsDict(db_dict).to_dict())

    else:
        msg = (
            f"unsupported output file extension {file_extension!r} for {args.output}: "
            "expected one of .json, .yaml, .yml, .pkl, .dat, .dir or .lh5"
        )
        raise ValueError(msg)
