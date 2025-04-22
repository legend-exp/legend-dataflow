import argparse
import json
import time
from pathlib import Path

from dbetto.catalog import Props
from legendmeta import TextDB
from lgdo import lh5
from pygama.hit.build_hit import build_hit

from ...alias_table import alias_table
from ...log import build_log


def build_tier_hit() -> None:
    argparser = argparse.ArgumentParser()
    argparser.add_argument("--input")
    argparser.add_argument("--pars-file", nargs="*")

    argparser.add_argument("--configs", required=True)
    argparser.add_argument("--table-map", required=False, type=str)
    argparser.add_argument("--log")
    argparser.add_argument("--alias-table", help="Alias table", type=str, default=None)

    argparser.add_argument("--datatype", required=True)
    argparser.add_argument("--timestamp", required=True)
    argparser.add_argument("--tier", required=True)

    argparser.add_argument("--output")
    args = argparser.parse_args()

    table_map = json.loads(args.table_map) if args.table_map is not None else None

    df_config = (
        TextDB(args.configs, lazy=True)
        .on(args.timestamp, system=args.datatype)
        .snakemake_rules[f"tier_{args.tier}"]
    )
    log = build_log(df_config, args.log, fallback=__name__)
    log.info("initializing")

    settings_dict = df_config.options.get("settings", {})

    if isinstance(settings_dict, str):
        settings_dict = Props.read_from(settings_dict)

    # mapping channel -> hit config file
    chan_cfg_map = df_config.inputs.hit_config

    log.info("building the build_hit config")
    # if the mapping only contains one __default__ key, build the channel
    # list from the (processable) channel map and assign the default config
    if list(chan_cfg_map.keys()) == ["__default__"]:
        chan_cfg_map = {chan: chan_cfg_map.__default__ for chan in table_map}

    # now construct the dictionary of hit configs for build_hit()
    channel_dict = {}
    pars_dict = {ch: chd["pars"] for ch, chd in Props.read_from(args.pars_file).items()}
    for chan, file in chan_cfg_map.items():
        hit_cfg = Props.read_from(file)

        # get pars (to override hit config)
        Props.add_to(hit_cfg, pars_dict.get(chan, {}).copy())

        if chan in table_map:
            input_tbl_name = table_map[chan] if table_map is not None else chan + "/dsp"
        else:
            continue

        # check if the raw tables are all existing
        if len(lh5.ls(args.input, input_tbl_name)) > 0:
            channel_dict[input_tbl_name] = hit_cfg
        else:
            msg = f"table {input_tbl_name} not found in {args.input} skipping"
            log.warning(msg)

    log.info("running build_hit()...")
    start = time.time()
    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    build_hit(args.input, lh5_tables_config=channel_dict, outfile=args.output)
    log.info(f"Hit built in {time.time() - start:.2f} seconds")
    if args.alias_table is not None:
        log.info("Creating alias table")
        alias_table(args.output, args.alias_table)


def build_tier_hit_single_channel() -> None:
    argparser = argparse.ArgumentParser()
    argparser.add_argument("--input")
    argparser.add_argument("--pars-file", nargs="*")

    argparser.add_argument("--configs", required=True)
    argparser.add_argument("--log")

    argparser.add_argument(
        "--channel",
        help="channel to process",
        required=False,
        type=str,
    )
    argparser.add_argument("--datatype", required=True)
    argparser.add_argument("--timestamp", required=True)
    argparser.add_argument("--tier", required=True)

    argparser.add_argument("--output")
    args = argparser.parse_args()

    df_config = (
        TextDB(args.configs, lazy=True)
        .on(args.timestamp, system=args.datatype)
        .snakemake_rules[f"tier_{args.tier}"]
    )
    log = build_log(df_config, args.log, fallback=__name__)
    log.info("initializing")

    settings_dict = df_config.options.get("settings", {})

    if isinstance(settings_dict, str):
        settings_dict = Props.read_from(settings_dict)

    # mapping channel -> hit config file
    chan_cfg_map = df_config.inputs.hit_config
    chan_cfg_map = (
        chan_cfg_map[args.channel]
        if args.channel is not None and args.channel in chan_cfg_map
        else chan_cfg_map
    )

    # now construct the dictionary of hit configs for build_hit()
    channel_dict = {}
    pars_dict = Props.read_from(args.pars_file)
    pars_dict = (
        pars_dict[args.channel]
        if args.channel is not None and args.channel in pars_dict
        else pars_dict
    )

    hit_cfg = Props.read_from(chan_cfg_map)
    Props.add_to(hit_cfg, pars_dict.copy())

    log.info("running build_hit()...")
    start = time.time()
    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    build_hit(args.input, hit_config=channel_dict, outfile=args.output)
    log.info(f"Hit built in {time.time() - start:.2f} seconds")
