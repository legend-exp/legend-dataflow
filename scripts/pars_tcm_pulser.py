import argparse
import logging
from pathlib import Path

import lgdo.lh5 as lh5
import numpy as np
from legendmeta import LegendMetadata, TextDB
from legendmeta.catalog import Props
from pygama.pargen.data_cleaning import get_tcm_pulser_ids

argparser = argparse.ArgumentParser()
argparser.add_argument("--configs", help="configs path", type=str, required=True)
argparser.add_argument("--metadata", help="metadata", type=str, required=True)
argparser.add_argument("--log", help="log file", type=str)

argparser.add_argument("--datatype", help="Datatype", type=str, required=True)
argparser.add_argument("--timestamp", help="Timestamp", type=str, required=True)
argparser.add_argument("--channel", help="Channel", type=str, required=True)

argparser.add_argument("--pulser_file", help="pulser file", type=str, required=False)

argparser.add_argument("--tcm_files", help="tcm_files", nargs="*", type=str)
args = argparser.parse_args()

configs = TextDB(args.configs, lazy=True).on(args.timestamp, system=args.datatype)
config_dict = configs["snakemake_rules"]["pars_tcm_pulser"]
if "logging" in config_dict["options"]:
    log_config = config_dict["options"]["logging"]
    log_config = Props.read_from(log_config)
    if args.log is not None:
        Path(args.log).parent.mkdir(parents=True, exist_ok=True)
        log_config["handlers"]["file"]["filename"] = args.log
    logging.config.dictConfig(log_config)
    log = logging.getLogger(config_dict["options"].get("logger", "prod"))
else:
    if args.log is not None:
        Path(args.log).parent.makedir(parents=True, exist_ok=True)
        logging.basicConfig(level=logging.INFO, filename=args.log, filemode="w")
    log = logging.getLogger(__name__)

sto = lh5.LH5Store()
log = logging.getLogger(__name__)


kwarg_dict = config_dict["inputs"]["pulser_config"]
kwarg_dict = Props.read_from(kwarg_dict)

meta = LegendMetadata(path=args.metadata)
channel_dict = meta.channelmap(args.timestamp, system=args.datatype)
channel = f"ch{channel_dict[args.channel].daq.rawid}"

if isinstance(args.tcm_files, list) and args.tcm_files[0].split(".")[-1] == "filelist":
    tcm_files = args.tcm_files[0]
    with Path(tcm_files).open() as f:
        tcm_files = f.read().splitlines()
else:
    tcm_files = args.tcm_files
# get pulser mask from tcm files
tcm_files = sorted(np.unique(tcm_files))
ids, mask = get_tcm_pulser_ids(tcm_files, channel, kwarg_dict.pop("pulser_multiplicity_threshold"))

Path(args.pulser_file).parent.mkdir(parents=True, exist_ok=True)
Props.write_to(args.pulser_file, {"idxs": ids.tolist(), "mask": mask.tolist()})
