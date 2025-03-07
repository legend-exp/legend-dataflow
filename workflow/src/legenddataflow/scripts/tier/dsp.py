import argparse
import re
from pathlib import Path

import numpy as np
from dbetto import TextDB
from dbetto.catalog import Props
from dspeed import build_dsp
from legendmeta import LegendMetadata
from lgdo import lh5

from ...log import build_log


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
    argparser.add_argument("--metadata", help="metadata repository path", required=True)
    argparser.add_argument("--log", help="log file name")

    argparser.add_argument("--datatype", help="datatype", required=True)
    argparser.add_argument("--timestamp", help="timestamp", required=True)
    argparser.add_argument("--tier", help="tier", required=True)

    argparser.add_argument(
        "--pars-file", help="database file for HPGes", nargs="*", default=[]
    )
    argparser.add_argument("--input", help="input file")

    argparser.add_argument("--output", help="output file")
    argparser.add_argument("--db-file", help="database file")
    args = argparser.parse_args()

    df_configs = TextDB(args.configs, lazy=True)
    config_dict = df_configs.on(args.timestamp, system=args.datatype).snakemake_rules

    if args.tier in ["dsp", "psp"]:
        config_dict = config_dict.tier_dsp
    elif args.tier in ["ann", "pan"]:
        config_dict = config_dict.tier_ann
    else:
        msg = f"tier {args.tier} not supported"
        raise ValueError(msg)

    log = build_log(config_dict, args.log, fallback=__name__)

    settings_dict = config_dict.options.get("settings", {})
    if isinstance(settings_dict, str):
        settings_dict = Props.read_from(settings_dict)

    chan_map = LegendMetadata(args.metadata).channelmap(
        args.timestamp, system=args.datatype
    )
    chan_cfg_map = config_dict.inputs.processing_chain

    # if the dictionary only contains one __default__ key, build the channel
    # list from the (processable) channel map and assign the default config
    if list(chan_cfg_map.keys()) == ["__default__"]:
        chan_cfg_map = {
            chan: chan_cfg_map.__default__
            for chan in chan_map.group("analysis.processable")[True].map("name")
        }

    # now construct the dictionary of DSP configs for build_dsp()
    dsp_cfg_tbl_dict = {}
    for chan, file in chan_cfg_map.items():
        if chan_map[chan].analysis.processable is False:
            msg = f"channel {chan} is set to non-processable in the channel map"
            raise RuntimeError(msg)

        tbl = "dsp" if args.tier in ["ann", "pan"] else "raw"
        input_tbl_name = f"ch{chan_map[chan].daq.rawid:07}/{tbl}"

        # check if the raw tables are all existing
        if len(lh5.ls(args.input, input_tbl_name)) > 0:
            dsp_cfg_tbl_dict[input_tbl_name] = Props.read_from(file)
        else:
            msg = f"table {input_tbl_name} not found in {args.input} skipping"
            log.warning(msg)

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
        (f"ch{chan_map[chan].daq.rawid:07}" if chan in chan_map else chan): dic
        for chan, dic in database_dict.items()
    }

    Path(args.output).parent.mkdir(parents=True, exist_ok=True)

    build_dsp(
        args.input,
        args.output,
        {},
        database=database_dict,
        chan_config=dsp_cfg_tbl_dict,
        write_mode="r",
        buffer_len=settings_dict.get("buffer_len", 1000),
        block_width=settings_dict.get("block_width", 16),
    )

    key = Path(args.output).name.replace(f"-tier_{args.tier}.lh5", "")

    if args.tier in ["dsp", "psp"]:
        raw_channels = [
            channel for channel in lh5.ls(args.input) if re.match("(ch\\d{7})", channel)
        ]
        raw_fields = [
            field.split("/")[-1]
            for field in lh5.ls(args.input, f"{raw_channels[0]}/raw/")
        ]

        outputs = {}
        channels = []
        for channel, chan_dict in dsp_cfg_tbl_dict.items():
            output = chan_dict["outputs"]
            in_dict = False
            for entry in outputs:
                if outputs[entry]["fields"] == output:
                    outputs[entry]["channels"].append(channel.split("/")[0])
                    in_dict = True
            if in_dict is False:
                outputs[f"group{len(list(outputs))+1}"] = {
                    "channels": [channel.split("/")[0]],
                    "fields": output,
                }
            channels.append(channel.split("/")[0])

        full_dict = {
            "valid_fields": {
                "raw": {"group1": {"fields": raw_fields, "channels": raw_channels}},
                "dsp": outputs,
            },
            "valid_keys": {
                key: {"valid_channels": {"raw": raw_channels, "dsp": channels}}
            },
        }
    else:
        outputs = {}
        channels = []
        for channel, chan_dict in dsp_cfg_tbl_dict.items():
            output = chan_dict["outputs"]
            in_dict = False
            for entry in outputs:
                if outputs[entry]["fields"] == output:
                    outputs[entry]["channels"].append(channel.split("/")[0])
                    in_dict = True
            if in_dict is False:
                outputs[f"group{len(list(outputs))+1}"] = {
                    "channels": [channel.split("/")[0]],
                    "fields": output,
                }
            channels.append(channel.split("/")[0])

        full_dict = {
            "valid_fields": {
                "ann": outputs,
            },
            "valid_keys": {key: {"valid_channels": {"ann": channels}}},
        }

    Path(args.db_file).parent.mkdir(parents=True, exist_ok=True)
    Props.write_to(args.db_file, full_dict)
