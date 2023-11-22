import argparse
import json
import logging
import os
import pathlib
import pickle as pkl
import time

import pygama.pargen.noise_optimization as pno
from legendmeta import LegendMetadata
from pygama.dsp.utils import numba_defaults

numba_defaults.cache = False
numba_defaults.boundscheck = True

argparser = argparse.ArgumentParser()
argparser.add_argument("--raw_filelist", help="raw_filelist", type=str)
argparser.add_argument("--database", help="database", type=str, required=True)
argparser.add_argument("--configs", help="configs", type=str, required=True)

argparser.add_argument("--log", help="log_file", type=str)

argparser.add_argument("--datatype", help="Datatype", type=str, required=True)
argparser.add_argument("--timestamp", help="Timestamp", type=str, required=True)
argparser.add_argument("--channel", help="Channel", type=str, required=True)

argparser.add_argument("--dsp_pars", help="dsp_pars", type=str, required=True)
argparser.add_argument("--plot_path", help="plot_path", type=str)

args = argparser.parse_args()

logging.basicConfig(level=logging.DEBUG, filename=args.log, filemode="w")
logging.getLogger("numba").setLevel(logging.INFO)
logging.getLogger("parse").setLevel(logging.INFO)
logging.getLogger("pygama.lgdo.lh5_store").setLevel(logging.INFO)
logging.getLogger("h5py._conv").setLevel(logging.INFO)
logging.getLogger("pygama.dsp.processing_chain").setLevel(logging.INFO)

log = logging.getLogger(__name__)


t0 = time.time()

conf = LegendMetadata(path=args.configs)
configs = conf.on(args.timestamp, system=args.datatype)
dsp_config = configs["snakemake_rules"]["pars_dsp_eopt"]["inputs"]["processing_chain"][
    args.channel
]
opt_json = configs["snakemake_rules"]["pars_dsp_nopt"]["inputs"]["optimiser_config"][args.channel]

with open(opt_json) as r:
    opt_dict = json.load(r)

with open(args.database) as t:
    db_dict = json.load(t)

if opt_dict["run_nopt"] is True:
    with open(args.raw_filelist) as f:
        files = f.read().splitlines()

    raw_files = sorted(files)

    if isinstance(dsp_config, str):
        with open(dsp_config) as r:
            dsp_config = json.load(r)

    if args.plot_path:
        out_dict, plot_dict = pno.noise_optimization(
            raw_files, dsp_config, db_dict, opt_dict, args.channel, display=1
        )
        pathlib.Path(os.path.dirname(args.plot_path)).mkdir(parents=True, exist_ok=True)
        with open(args.plot_path, "wb") as f:
            pkl.dump(plot_dict, f, protocol=pkl.HIGHEST_PROTOCOL)
    else:
        out_dict = pno.noise_optimization(raw_files, dsp_config, db_dict, opt_dict, args.channel)

    t1 = time.time()
    log.info(f"Optimiser finished in {(t1-t0)/60} minutes")
else:
    out_dict = {}

pathlib.Path(os.path.dirname(args.dsp_pars)).mkdir(parents=True, exist_ok=True)
with open(args.dsp_pars, "w") as w:
    json.dump(db_dict, w, indent=4)
