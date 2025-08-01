from __future__ import annotations

import argparse
import time
from pathlib import Path

import hdf5plugin
from daq2lh5 import build_raw
from dbetto import TextDB
from dbetto.catalog import Props
from legenddataflowscripts.utils import alias_table, build_log

filter_map = {
    "zstd": hdf5plugin.Zstd(),
    "blosc": hdf5plugin.Blosc(),
    "lz4": hdf5plugin.LZ4(),
    "bitshuffle": hdf5plugin.Bitshuffle(),
    "bzip2": hdf5plugin.BZip2(),
    # Add other filters from hdf5plugin if needed
}


def build_tier_raw_fcio() -> None:
    argparser = argparse.ArgumentParser()
    argparser.add_argument("input", help="input file", type=str)
    argparser.add_argument("output", help="output file", type=str)
    argparser.add_argument("--datatype", help="Datatype", type=str, required=True)
    argparser.add_argument("--timestamp", help="Timestamp", type=str, required=True)
    argparser.add_argument("--configs", help="config file", type=str)
    argparser.add_argument("--chan-maps", help="chan map", type=str)
    argparser.add_argument("--log", help="log file", type=str)
    argparser.add_argument("--alias-table", help="Alias table", type=str, default=None)
    argparser.add_argument(
        "--inl-table", help="INL table", type=str, default=None, required=False
    )
    args = argparser.parse_args()

    Path(args.output).parent.mkdir(parents=True, exist_ok=True)

    config_dict = (
        TextDB(args.configs, lazy=True)
        .on(args.timestamp, system=args.datatype)
        .snakemake_rules.tier_raw_fcio
    )

    log = build_log(config_dict, args.log)

    channel_dict = config_dict.inputs
    settings = Props.read_from(channel_dict.settings)
    channel_dict = channel_dict.out_spec
    all_config = Props.read_from(channel_dict.gen_config)

    if "hdf5_settings" in settings and "compression" in settings["hdf5_settings"]:
        compression = settings["hdf5_settings"]["compression"]
        if compression in filter_map:
            settings["hdf5_settings"]["compression"] = filter_map[compression]

    if "geds_config" in channel_dict:
        raise NotImplementedError()

    if "spms_config" in channel_dict:
        spm_config = Props.read_from(channel_dict.spms_config)
        Props.add_to(all_config, spm_config)

    if "auxs_config" in channel_dict:
        raise NotImplementedError()

    if "muon_config" in channel_dict:
        raise NotImplementedError()

    log.info("Built raw config")
    log.info("Starting build raw")
    start = time.time()

    build_raw(
        args.input,
        out_spec=all_config,
        in_stream_type="FlashCam",
        filekey=args.output,
        **settings,
    )

    msg = f"Built raw in {time.time() - start:.2f} seconds"
    log.info(msg)

    if args.alias_table is not None:
        log.info("Creating alias table")
        alias_table(args.output, args.alias_table)
