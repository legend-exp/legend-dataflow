import argparse
import logging
import os
import pathlib

os.environ["LGDO_CACHE"] = "false"
os.environ["LGDO_BOUNDSCHECK"] = "false"
os.environ["PYGAMA_PARALLEL"] = "false"
os.environ["PYGAMA_FASTMATH"] = "false"

import lgdo.lh5 as lh5
import numpy as np
from legendmeta import LegendMetadata
from legendmeta.catalog import Props
from pygama.pargen.data_cleaning import get_tcm_pulser_ids

argparser = argparse.ArgumentParser()
argparser.add_argument("--configs", help="configs path", type=str, required=True)
argparser.add_argument("--log", help="log file", type=str)

argparser.add_argument("--datatype", help="Datatype", type=str, required=True)
argparser.add_argument("--timestamp", help="Timestamp", type=str, required=True)
argparser.add_argument("--channel", help="Channel", type=str, required=True)

argparser.add_argument("--pulser_file", help="pulser file", type=str, required=False)

argparser.add_argument("--tcm_files", help="tcm_files", nargs="*", type=str)
args = argparser.parse_args()

logging.basicConfig(level=logging.DEBUG, filename=args.log, filemode="w")
logging.getLogger("numba").setLevel(logging.INFO)
logging.getLogger("parse").setLevel(logging.INFO)
logging.getLogger("lgdo").setLevel(logging.INFO)
logging.getLogger("h5py").setLevel(logging.INFO)
logging.getLogger("matplotlib").setLevel(logging.INFO)
logging.getLogger("legendmeta").setLevel(logging.INFO)

sto = lh5.LH5Store()
log = logging.getLogger(__name__)

configs = LegendMetadata(path=args.configs)
config_dict = configs.on(args.timestamp, system=args.datatype)
kwarg_dict = config_dict["snakemake_rules"]["pars_tcm_pulser"]["inputs"]["pulser_config"]

kwarg_dict = Props.read_from(kwarg_dict)

if isinstance(args.tcm_files, list) and args.tcm_files[0].split(".")[-1] == "filelist":
    tcm_files = args.tcm_files[0]
    with open(tcm_files) as f:
        tcm_files = f.read().splitlines()
else:
    tcm_files = args.tcm_files
# get pulser mask from tcm files
tcm_files = sorted(np.unique(tcm_files))
ids, mask = get_tcm_pulser_ids(
    tcm_files, args.channel, kwarg_dict.pop("pulser_multiplicity_threshold")
)

pathlib.Path(os.path.dirname(args.pulser_file)).mkdir(parents=True, exist_ok=True)
Props.write_to(args.pulser_file, {"idxs": ids.tolist(), "mask": mask.tolist()})
