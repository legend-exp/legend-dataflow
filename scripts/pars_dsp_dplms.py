from pygama.dsp.utils import numba_defaults

numba_defaults.cache = False
numba_defaults.boundscheck = True

import argparse
import json
import logging
import os
import pathlib
import pickle as pkl
import time

import pygama.pargen.dplms_ge_dict as pdd
from legendmeta import LegendMetadata

argparser = argparse.ArgumentParser()
argparser.add_argument("--fft_raw_filelist", help="fft_raw_filelist", type=str)
argparser.add_argument("--cal_raw_filelist", help="cal_raw_filelist", type=str)
argparser.add_argument("--database", help="database", type=str, required=True)
argparser.add_argument("--configs", help="configs", type=str, required=True)

argparser.add_argument("--log", help="log_file", type=str)

argparser.add_argument("--datatype", help="Datatype", type=str, required=True)
argparser.add_argument("--timestamp", help="Timestamp", type=str, required=True)
argparser.add_argument("--channel", help="Channel", type=str, required=True)

argparser.add_argument("--dsp_pars", help="dsp_pars", type=str, required=True)
argparser.add_argument("--lh5_path", help="lh5_path", type=str, required=True)
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
configs = configs.on(args.timestamp, system=args.datatype)
dsp_config = config_dict['snakemake_rules']['pars_dsp_dplms']["inputs"]['proc_chain'][args.channel]

dplms_json = config_dict['snakemake_rules']['pars_dsp_dplms']["inputs"]['dplms_pars'][args.channel]
with open(dplms_json) as r:
    dplms_dict = json.load(r)

with open(args.database) as t:
    db_dict = json.load(t)

if opt_dict["run_dplms"] is True:
    with open(args.fft_raw_filelist) as f:
        fft_files = f.read().splitlines()
    with open(args.cal_raw_filelist) as f:
        cal_files = f.read().splitlines()

    fft_files = sorted(fft_files)
    cal_files = sorted(cal_files)

    if isinstance(dsp_config, str):
        with open(dsp_config) as r:
            dsp_config = json.load(r)

    if args.plot_path:
        out_dict, plot_dict = pdd.dplms_ge_dict(
            args.channel,
            fft_files,
            cal_files,
            dsp_config,
            db_dict,
            args.lh5_path,
            dplms_dict,
            display=1
        )
        pathlib.Path(os.path.dirname(args.plot_path)).mkdir(parents=True, exist_ok=True)
        with open(args.plot_path, "wb") as f:
            pkl.dump(plot_dict, f, protocol=pkl.HIGHEST_PROTOCOL)
    else:
        out_dict, plot_dict = pdd.dplms_ge_dict(
            args.channel,
            fft_files,
            cal_files,
            dsp_config,
            db_dict,
            args.lh5_path,
            dplms_dict,
        )

    t1 = time.time()
    log.info(f"DPLMS creation finished in {(t1-t0)/60} minutes")
else:
    out_dict = {}

pathlib.Path(os.path.dirname(args.dsp_pars)).mkdir(parents=True, exist_ok=True)
with open(args.dsp_pars, "w") as w:
    json.dump(out_dict, w, indent=2)
