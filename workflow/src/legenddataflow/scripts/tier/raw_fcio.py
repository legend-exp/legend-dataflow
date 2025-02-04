import argparse
from copy import deepcopy
from pathlib import Path

import numpy as np
from daq2lh5 import build_raw
from dbetto import TextDB
from dbetto.catalog import Props

from ...log import build_log


def build_tier_raw_fcio() -> None:
    argparser = argparse.ArgumentParser()
    argparser.add_argument("input", help="input file", type=str)
    argparser.add_argument("output", help="output file", type=str)
    argparser.add_argument("--datatype", help="Datatype", type=str, required=True)
    argparser.add_argument("--timestamp", help="Timestamp", type=str, required=True)
    argparser.add_argument("--configs", help="config file", type=str)
    argparser.add_argument("--chan-maps", help="chan map", type=str)
    argparser.add_argument("--log", help="log file", type=str)
    args = argparser.parse_args()

    Path(args.output).parent.mkdir(parents=True, exist_ok=True)

    config_dict = (
        TextDB(args.configs, lazy=True)
        .on(args.timestamp, system=args.datatype)
        .snakemake_rules.tier_raw_fcio
    )

    build_log(config_dict, args.log)

    channel_dict = config_dict.inputs
    settings = Props.read_from(channel_dict.settings)
    channel_dict = channel_dict.out_spec
    all_config = Props.read_from(channel_dict.gen_config)

    chmap = (
        TextDB(args.chan_maps, lazy=True).channelmaps.on(args.timestamp).group("system")
    )

    if "geds_config" in channel_dict:
        raise NotImplementedError()

    if "spms_config" in channel_dict:
        spm_config = Props.read_from(channel_dict.spms_config)
        spm_channels = chmap.spms.map("daq.rawid")

        for rawid, chinfo in spm_channels.items():
            cfg_block = deepcopy(spm_config["FCEventDecoder"]["__output_table_name__"])
            cfg_block["key_list"] = [chinfo.daq.fc_channel]
            spm_config["FCEventDecoder"][f"ch{rawid:07d}/raw"] = cfg_block

        spm_config["FCEventDecoder"].pop("__output_table_name__")

        Props.add_to(all_config, spm_config)

    if "auxs_config" in channel_dict:
        raise NotImplementedError()

    if "muon_config" in channel_dict:
        raise NotImplementedError()

    rng = np.random.default_rng()
    rand_num = f"{rng.integers(0,99999):05d}"
    temp_output = f"{args.output}.{rand_num}"

    build_raw(args.input, out_spec=all_config, filekey=temp_output, **settings)

    # rename the temp file
    Path(temp_output).rename(args.output)
