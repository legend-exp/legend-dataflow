import argparse
import logging
import time
from pathlib import Path

from daq2lh5 import build_raw
from dbetto import TextDB
from dbetto.catalog import Props

from ...alias_table import alias_table
from ...log import build_log


def build_tier_raw_orca() -> None:
    argparser = argparse.ArgumentParser()
    argparser.add_argument("input", help="input file", type=str)
    argparser.add_argument("output", help="output file", type=str)
    argparser.add_argument("--datatype", help="Datatype", type=str, required=True)
    argparser.add_argument("--timestamp", help="Timestamp", type=str, required=True)
    argparser.add_argument("--configs", help="config file", type=str)
    argparser.add_argument("--chan-maps", help="chan map", type=str)
    argparser.add_argument("--log", help="log file")
    argparser.add_argument("--alias-table", help="Alias table", type=str, default=None)
    args = argparser.parse_args()

    Path(args.log).parent.mkdir(parents=True, exist_ok=True)
    logging.basicConfig(level=logging.INFO, filename=args.log, filemode="w")

    Path(args.output).parent.mkdir(parents=True, exist_ok=True)

    configs = TextDB(args.configs, lazy=True)
    config_dict = configs.on(args.timestamp, system=args.datatype)["snakemake_rules"][
        "tier_raw_orca"
    ]

    log = build_log(config_dict, args.log)

    channel_dict = config_dict["inputs"]
    settings = Props.read_from(channel_dict["settings"])
    channel_dict = channel_dict["out_spec"]
    all_config = Props.read_from(channel_dict["gen_config"])

    chmap = TextDB(args.chan_maps, lazy=True).channelmaps.on(args.timestamp)

    if "geds_config" in list(channel_dict):
        ged_config = Props.read_from(channel_dict["geds_config"])

        ged_channels = list(chmap.map("system", unique=False)["geds"].map("daq.rawid"))

        ged_config[next(iter(ged_config))]["geds"]["key_list"] = sorted(ged_channels)
        Props.add_to(all_config, ged_config)

    if "spms_config" in list(channel_dict):
        spm_config = Props.read_from(channel_dict["spms_config"])

        spm_channels = list(chmap.map("system", unique=False)["spms"].map("daq.rawid"))

        spm_config[next(iter(spm_config))]["spms"]["key_list"] = sorted(spm_channels)
        Props.add_to(all_config, spm_config)

    if "muon_config" in list(channel_dict):
        muon_config = Props.read_from(channel_dict["muon_config"])
        muon_channels = list(chmap.map("system", unique=False)["muon"].map("daq.rawid"))
        top_key = next(iter(muon_config))
        muon_config[top_key][next(iter(muon_config[top_key]))]["key_list"] = sorted(
            muon_channels
        )
        Props.add_to(all_config, muon_config)

    if "auxs_config" in list(channel_dict):
        aux_config = Props.read_from(channel_dict["auxs_config"])
        aux_channels = list(
            chmap.map("system", unique=False)["auxs"]
            .map("daq.crate", unique=False)[1]
            .map("daq.rawid")
        )
        aux_channels += list(
            chmap.map("system", unique=False)["auxs"]
            .map("daq.crate", unique=False)[0]
            .map("daq.rawid")
        )
        aux_channels += list(chmap.map("system", unique=False)["puls"].map("daq.rawid"))
        aux_channels += list(chmap.map("system", unique=False)["bsln"].map("daq.rawid"))
        top_key = next(iter(aux_config))
        aux_config[top_key][next(iter(aux_config[top_key]))]["key_list"] = sorted(
            aux_channels
        )
        Props.add_to(all_config, aux_config)

    log.info("Built raw config")
    log.info("Starting build raw")
    start = time.time()

    build_raw(
        args.input,
        out_spec=all_config,
        in_stream_type="ORCA",
        filekey=args.output,
        **settings,
    )

    log.info(f"Built raw in {time.time() - start:.2f} seconds")

    if args.alias_table is not None:
        log.info("Creating alias table")
        alias_table(args.output, args.alias_table)
