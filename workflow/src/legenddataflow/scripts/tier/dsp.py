import argparse
import json
import time
import warnings
from pathlib import Path

import numpy as np
from dbetto import TextDB
from dbetto.catalog import Props
from dspeed import build_dsp
from lgdo import lh5

from ...alias_table import alias_table
from ...log import build_log

warnings.filterwarnings(action="ignore", category=RuntimeWarning)


def _replace_list_with_array(dic):
    for key, value in dic.items():
        if isinstance(value, dict):
            dic[key] = _replace_list_with_array(value)
        elif isinstance(value, list):
            dic[key] = np.array(value, dtype="float32")
        else:
            pass
    return dic


def build_tier_dsp() -> None:
    # CLI config
    argparser = argparse.ArgumentParser()
    argparser.add_argument(
        "--configs", help="path to dataflow config files", required=True
    )
    argparser.add_argument(
        "--table-map",
        help="mapping from channel to table name",
        required=False,
        type=str,
    )
    argparser.add_argument("--log", help="log file name")
    argparser.add_argument("--alias-table", help="Alias table", type=str, default=None)

    argparser.add_argument("--n-processes", help="log file name", default=1, type=int)

    argparser.add_argument("--datatype", help="datatype", required=True)
    argparser.add_argument("--timestamp", help="timestamp", required=True)
    argparser.add_argument("--tier", help="tier", required=True)

    argparser.add_argument(
        "--pars-file", help="database file for HPGes", nargs="*", default=[]
    )
    argparser.add_argument("--input", help="input file")

    argparser.add_argument("--output", help="output file")
    args = argparser.parse_args()

    # set number of threads to use
    # set_num_threads(1)

    table_map = json.loads(args.table_map) if args.table_map is not None else None

    df_configs = TextDB(args.configs, lazy=True)
    config_dict = df_configs.on(args.timestamp, system=args.datatype).snakemake_rules
    config_dict = config_dict[f"tier_{args.tier}"]

    log = build_log(config_dict, args.log, fallback=__name__)

    settings_dict = config_dict.options.get("settings", {})
    if isinstance(settings_dict, str):
        settings_dict = Props.read_from(settings_dict)

    chan_cfg_map = config_dict.inputs.processing_chain

    # if the dictionary only contains one __default__ key, build the channel
    # list from the (processable) channel map and assign the default config
    if list(chan_cfg_map.keys()) == ["__default__"]:
        chan_cfg_map = {chan: chan_cfg_map.__default__ for chan in table_map}

    # now construct the dictionary of DSP configs for build_dsp()
    dsp_cfg_tbl_dict = {}
    for chan, file in chan_cfg_map.items():
        if chan in table_map:
            input_tbl_name = table_map[chan] if table_map is not None else chan + "/raw"
        else:
            continue

        # check if the raw tables are all existing
        if len(lh5.ls(args.input, input_tbl_name)) > 0:
            dsp_cfg_tbl_dict[input_tbl_name] = Props.read_from(file)
        else:
            msg = f"table {input_tbl_name} not found in {args.input} skipping"
            log.info(msg)

    if len(dsp_cfg_tbl_dict) == 0:
        msg = f"could not find any of the requested channels in {args.input}"
        raise RuntimeError(msg)

    # par files
    db_files = [
        par_file
        for par_file in args.pars_file
        if Path(par_file).suffix in (".json", ".yaml", ".yml")
    ]

    database_dict = _replace_list_with_array(
        Props.read_from(db_files, subst_pathvar=True)
    )
    database_dict = {
        (table_map[chan].split("/")[0] if chan in table_map else chan): dic
        for chan, dic in database_dict.items()
    }
    log.info("loaded database files")

    Path(args.output).parent.mkdir(parents=True, exist_ok=True)

    start = time.time()

    build_dsp(
        args.input,
        args.output,
        database=database_dict,
        chan_config=dsp_cfg_tbl_dict,
        write_mode="r",
        buffer_len=settings_dict.get("buffer_len", 1000),
        block_width=settings_dict.get("block_width", 16),
        # n_processes=args.n_processes,
    )

    log.info(f"Finished building DSP in {time.time()- start:.2f} seconds")
    if args.alias_table is not None:
        log.info("Creating alias table")
        alias_table(args.output, args.alias_table)


def build_tier_dsp_single_channel() -> None:
    # CLI config
    argparser = argparse.ArgumentParser()
    argparser.add_argument(
        "--configs", help="path to dataflow config files", required=True
    )
    argparser.add_argument(
        "--channel",
        help="channel to process",
        required=False,
        type=str,
    )
    argparser.add_argument("--log", help="log file name")

    argparser.add_argument("--datatype", help="datatype", required=True)
    argparser.add_argument("--timestamp", help="timestamp", required=True)
    argparser.add_argument("--tier", help="tier", required=True)

    argparser.add_argument(
        "--pars-file", help="database file for HPGes", nargs="*", default=[]
    )
    argparser.add_argument("--input", help="input file")

    argparser.add_argument("--output", help="output file")
    args = argparser.parse_args()

    df_configs = TextDB(args.configs, lazy=True)
    config_dict = df_configs.on(args.timestamp, system=args.datatype).snakemake_rules
    config_dict = config_dict[f"tier_{args.tier}"]
    config_dict = (
        config_dict[args.channel]
        if args.channel is not None and args.channel in config_dict
        else config_dict
    )

    log = build_log(config_dict, args.log, fallback=__name__)

    settings_dict = config_dict.options.get("settings", {})
    if isinstance(settings_dict, str):
        settings_dict = Props.read_from(settings_dict)

    proc_chain = config_dict.inputs.processing_chain

    # par files
    db_files = [
        par_file
        for par_file in args.pars_file
        if Path(par_file).suffix in (".json", ".yaml", ".yml")
    ]

    database_dict = _replace_list_with_array(
        Props.read_from(db_files, subst_pathvar=True)
    )
    database_dict = (
        database_dict[args.channel]
        if args.channel is not None and args.channel in database_dict
        else database_dict
    )

    Path(args.output).parent.mkdir(parents=True, exist_ok=True)

    start = time.time()

    build_dsp(
        args.input,
        args.output,
        proc_chain,
        database=database_dict,
        write_mode="r",
        buffer_len=settings_dict.get("buffer_len", 1000),
        block_width=settings_dict.get("block_width", 16),
    )
    log.info(f"Finished building DSP in {time.time()- start:.2f} seconds")
