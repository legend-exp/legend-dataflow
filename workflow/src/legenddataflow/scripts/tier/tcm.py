import argparse

import lgdo.lh5 as lh5
from daq2lh5.orca import orca_flashcam
from dbetto import TextDB
from dbetto.catalog import Props
from pygama.evt.build_tcm import build_tcm

from ...log import build_log


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

    build_log(config_dict, args.log)

    settings = Props.read_from(config_dict["inputs"]["config"])

    # get the list of channels by fcid
    ch_list = lh5.ls(args.input, "/ch*")
    fcid_channels = {}
    for ch in ch_list:
        key = int(ch[2:])
        fcid = orca_flashcam.get_fcid(key)
        if fcid not in fcid_channels:
            fcid_channels[fcid] = []
        fcid_channels[fcid].append(f"/{ch}/raw")

    # make a hardware_tcm_[fcid] for each fcid
    for fcid, fcid_dict in fcid_channels.items():
        build_tcm(
            [(args.input, fcid_dict)],
            out_file=args.output,
            out_name=f"hardware_tcm_{fcid}",
            wo_mode="o",
            **settings,
        )
