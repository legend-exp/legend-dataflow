"""Console script ``build-tier-tcm``: build the time-coincidence map (tcm)
tier from raw-tier data."""

from __future__ import annotations

import argparse
from pathlib import Path

import lh5
from daq2lh5.orca import orca_flashcam
from dbetto import AttrsDict, TextDB
from dbetto.catalog import Props
from legenddataflowscripts.utils import build_log
from pygama.evt.build_tcm import build_tcm


def build_tier_tcm() -> None:
    argparser = argparse.ArgumentParser()
    argparser.add_argument("input", help="input file", type=str)
    argparser.add_argument("output", help="output file", type=str)
    argparser.add_argument("--datatype", help="Datatype", type=str, required=True)
    argparser.add_argument("--timestamp", help="Timestamp", type=str, required=True)
    argparser.add_argument("--configs", help="config file", type=str)
    argparser.add_argument("--log", help="log file", type=str)
    args = argparser.parse_args()

    configs = TextDB(args.configs, lazy=True).on(args.timestamp, system=args.datatype)
    config_dict = configs["snakemake_rules"]["tier_tcm"]

    log = build_log(config_dict, args.log)
    # the config can either be a single file used for all fcids or a mapping
    # of per-fcid config files
    per_fcid_settings = None
    if isinstance(config_dict["inputs"]["config"], dict | AttrsDict):
        per_fcid_settings = {
            int(key): Props.read_from(val)
            for key, val in config_dict["inputs"]["config"].items()
        }
        settings = None
    else:
        settings = Props.read_from(config_dict["inputs"]["config"])

    # get the list of channels by fcid
    ch_list = lh5.ls(args.input, "/ch*")

    if len(ch_list) == 0:
        msg = f"no tables matching /ch* found in input file {args.input}"
        raise RuntimeError(msg)

    log.debug(ch_list)

    fcid_channels = {}
    for ch in ch_list:
        key = int(ch[2:])
        fcid = orca_flashcam.get_fcid(key)
        if fcid not in fcid_channels:
            fcid_channels[fcid] = []
        fcid_channels[fcid].append(f"/{ch}/raw")

    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    # make a hardware_tcm_[fcid] for each fcid
    for fcid, fcid_dict in fcid_channels.items():
        msg = f"building tcm for fcid: {fcid}"
        log.info(msg)
        if per_fcid_settings is not None:
            if fcid not in per_fcid_settings:
                msg = (
                    f"no tcm config found for fcid {fcid}: the config mapping "
                    f"only covers fcids {sorted(per_fcid_settings)}"
                )
                raise RuntimeError(msg)
            fcid_settings = per_fcid_settings[fcid]
        else:
            fcid_settings = settings
        build_tcm(
            [(args.input, fcid_dict)],
            out_file=args.output,
            out_name=f"hardware_tcm_{fcid}",
            wo_mode="o",
            **fcid_settings,
        )
        msg = f"built tcm for fcid: {fcid}"
        log.info(msg)
