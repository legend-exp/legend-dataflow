import argparse
from pathlib import Path

import numpy as np
from dbetto import TextDB
from dbetto.catalog import Props
from legendmeta import LegendMetadata
from pygama.pargen.data_cleaning import get_tcm_pulser_ids

from ....log import build_log


def par_geds_tcm_pulser() -> None:
    argparser = argparse.ArgumentParser()
    argparser.add_argument("--configs", help="configs path", type=str, required=True)
    argparser.add_argument("--metadata", help="metadata", type=str, required=True)
    argparser.add_argument("--log", help="log file", type=str)

    argparser.add_argument("--datatype", help="Datatype", type=str, required=True)
    argparser.add_argument("--timestamp", help="Timestamp", type=str, required=True)
    argparser.add_argument("--channel", help="Channel", type=str, required=True)

    argparser.add_argument(
        "--pulser-file", help="pulser file", type=str, required=False
    )

    argparser.add_argument("--tcm-files", help="tcm_files", nargs="*", type=str)
    args = argparser.parse_args()

    configs = TextDB(args.configs, lazy=True).on(args.timestamp, system=args.datatype)
    config_dict = configs["snakemake_rules"]["pars_tcm_pulser"]

    build_log(config_dict, args.log)

    kwarg_dict = config_dict["inputs"]["pulser_config"]
    kwarg_dict = Props.read_from(kwarg_dict)

    meta = LegendMetadata(path=args.metadata)
    channel_dict = meta.channelmap(args.timestamp, system=args.datatype)
    channel = f"ch{channel_dict[args.channel].daq.rawid}"

    if (
        isinstance(args.tcm_files, list)
        and args.tcm_files[0].split(".")[-1] == "filelist"
    ):
        tcm_files = args.tcm_files[0]
        with Path(tcm_files).open() as f:
            tcm_files = f.read().splitlines()
    else:
        tcm_files = args.tcm_files
    # get pulser mask from tcm files
    tcm_files = sorted(np.unique(tcm_files))
    ids, mask = get_tcm_pulser_ids(
        tcm_files, channel, kwarg_dict.pop("pulser_multiplicity_threshold")
    )

    Path(args.pulser_file).parent.mkdir(parents=True, exist_ok=True)
    Props.write_to(args.pulser_file, {"idxs": ids.tolist(), "mask": mask.tolist()})
