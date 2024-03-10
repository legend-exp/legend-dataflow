import argparse
import json
import logging
import os
import pathlib
import pickle as pkl
import time

os.environ["LGDO_CACHE"] = "false"
os.environ["LGDO_BOUNDSCHECK"] = "false"
os.environ["DSPEED_CACHE"] = "false"
os.environ["DSPEED_BOUNDSCHECK"] = "false"

import lgdo.lh5 as lh5
import numpy as np
import pygama.pargen.noise_optimization as pno
from legendmeta import LegendMetadata
from legendmeta.catalog import Props
from pygama.pargen.cuts import generate_cuts, get_cut_indexes
from pygama.pargen.dsp_optimize import run_one_dsp

sto = lh5.LH5Store()

argparser = argparse.ArgumentParser()
argparser.add_argument("--raw_filelist", help="raw_filelist", type=str)
argparser.add_argument("--database", help="database", type=str, required=True)
argparser.add_argument("--inplots", help="inplots", type=str)

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
logging.getLogger("lgdo").setLevel(logging.INFO)
logging.getLogger("h5py._conv").setLevel(logging.INFO)
logging.getLogger("dspeed.processing_chain").setLevel(logging.INFO)
logging.getLogger("legendmeta").setLevel(logging.INFO)

log = logging.getLogger(__name__)


t0 = time.time()

conf = LegendMetadata(path=args.configs)
configs = conf.on(args.timestamp, system=args.datatype)
dsp_config = configs["snakemake_rules"]["pars_dsp_nopt"]["inputs"]["processing_chain"][
    args.channel
]
opt_json = configs["snakemake_rules"]["pars_dsp_nopt"]["inputs"]["optimiser_config"][args.channel]

opt_dict = Props.read_from(opt_json)

db_dict = Props.read_from(args.database)

if opt_dict.pop("run_nopt") is True:
    with open(args.raw_filelist) as f:
        files = f.read().splitlines()

    raw_files = sorted(files)

    energies = sto.read(f"{args.channel}/raw/daqenergy", raw_files)[0]
    idxs = np.where(energies.nda == 0)[0]
    tb_data = sto.read(f"{args.channel}/raw", raw_files, n_rows=opt_dict["n_events"], idx=idxs)[0]
    t1 = time.time()
    log.info(f"Time to open raw files {t1-t0:.2f} s, n. baselines {len(tb_data)}")

    log.info(f"Select baselines {len(tb_data)}")
    dsp_data = run_one_dsp(tb_data, dsp_config)
    cut_dict = generate_cuts(dsp_data, parameters=opt_dict.pop("cut_pars"))
    cut_idxs = get_cut_indexes(dsp_data, cut_dict)
    tb_data = sto.read(
        f"{args.channel}/raw", raw_files, n_rows=opt_dict.pop("n_events"), idx=idxs[cut_idxs]
    )[0]
    log.info(f"... {len(tb_data)} baselines after cuts")

    if isinstance(dsp_config, (str, list)):
        dsp_config = Props.read_from(dsp_config)

    if args.plot_path:
        out_dict, plot_dict = pno.noise_optimization(
            tb_data, dsp_config, db_dict.copy(), opt_dict, args.channel, display=1
        )
    else:
        out_dict = pno.noise_optimization(
            raw_files, dsp_config, db_dict.copy(), opt_dict, args.channel
        )

    t2 = time.time()
    log.info(f"Optimiser finished in {(t2-t0)/60} minutes")
else:
    out_dict = {}
    plot_dict = {}

if args.plot_path:
    pathlib.Path(os.path.dirname(args.plot_path)).mkdir(parents=True, exist_ok=True)
    if args.inplots:
        with open(args.inplots, "rb") as r:
            old_plot_dict = pkl.load(r)
        plot_dict = dict(noise_optimisation=plot_dict, **old_plot_dict)
    else:
        plot_dict = {"noise_optimisation": plot_dict}
    with open(args.plot_path, "wb") as f:
        pkl.dump(plot_dict, f, protocol=pkl.HIGHEST_PROTOCOL)

pathlib.Path(os.path.dirname(args.dsp_pars)).mkdir(parents=True, exist_ok=True)
with open(args.dsp_pars, "w") as w:
    json.dump(dict(nopt_pars=out_dict, **db_dict), w, indent=4)
