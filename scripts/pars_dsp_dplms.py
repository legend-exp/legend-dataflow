import argparse
import json
import logging
import os
import pathlib
import pickle as pkl
import time

import lgdo.lh5_store as lh5
import numpy as np
import pygama.pargen.dplms_ge_dict as pdd
from legendmeta import LegendMetadata
from pygama.dsp.utils import numba_defaults
from pygama.pargen.energy_optimisation import (
    event_selection,
    index_data,
)

numba_defaults.cache = False
numba_defaults.boundscheck = True

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
sto = lh5.LH5Store()

conf = LegendMetadata(path=args.configs)
configs = conf.on(args.timestamp, system=args.datatype)
dsp_config = configs["snakemake_rules"]["pars_dsp_dplms"]["inputs"]["proc_chain"][args.channel]

dplms_json = configs["snakemake_rules"]["pars_dsp_dplms"]["inputs"]["dplms_pars"][args.channel]
with open(dplms_json) as r:
    dplms_dict = json.load(r)

with open(args.database) as t:
    db_dict = json.load(t)

if dplms_dict["run_dplms"] is True:
    with open(args.fft_raw_filelist) as f:
        fft_files = f.read().splitlines()
    with open(args.cal_raw_filelist) as f:
        cal_files = f.read().splitlines()

    fft_files = sorted(fft_files)
    cal_files = sorted(cal_files)

    t0 = time.time()
    log.info("\nLoad fft data")
    energies = sto.read_object(f"{args.channel}/raw/daqenergy", fft_files)[0]
    idxs = np.where(energies.nda == 0)[0]
    raw_fft = sto.read_object(
        f"{args.channel}/raw", fft_files, n_rows=dplms_dict["n_baselines"], idx=idxs
    )[0]
    t1 = time.time()
    log.info(f"Time to load fft data {(t1-t0):.2f} s, total events {len(raw_fft)}")

    log.info("\nRunning event selection")
    peaks_keV = np.array(dplms_dict["peaks_keV"])
    kev_widths = [tuple(kev_width) for kev_width in dplms_dict["kev_widths"]]
    raw_cal, idx_list = event_selection(
        cal_files,
        f"{args.channel}/raw",
        dsp_config,
        db_dict[args.channel],
        peaks_keV,
        np.arange(0, len(peaks_keV), 1).tolist(),
        kev_widths,
        cut_parameters=dplms_dict["wfs_cut_pars"],
        n_events=dplms_dict["n_signals"],
    )
    raw_cal = index_data(raw_cal, idx_list[-1])
    log.info(f"Time to run event selection {(time.time()-t1):.2f} s, total events {len(raw_cal)}")

    if isinstance(dsp_config, str):
        with open(dsp_config) as r:
            dsp_config = json.load(r)

    if args.plot_path:
        out_dict, plot_dict = pdd.dplms_ge_dict(
            args.channel,
            raw_fft,
            raw_cal,
            dsp_config,
            db_dict,
            args.lh5_path,
            dplms_dict,
            display=1,
        )
        pathlib.Path(os.path.dirname(args.plot_path)).mkdir(parents=True, exist_ok=True)
        with open(args.plot_path, "wb") as f:
            pkl.dump(plot_dict, f, protocol=pkl.HIGHEST_PROTOCOL)
    else:
        out_dict, plot_dict = pdd.dplms_ge_dict(
            args.channel,
            raw_fft,
            raw_cal,
            dsp_config,
            db_dict,
            args.lh5_path,
            dplms_dict,
        )

    log.info(f"DPLMS creation finished in {(time.time()-t0)/60} minutes")
else:
    out_dict = {}

pathlib.Path(os.path.dirname(args.dsp_pars)).mkdir(parents=True, exist_ok=True)
with open(args.dsp_pars, "w") as w:
    json.dump(out_dict, w, indent=2)
