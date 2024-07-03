import argparse
import logging
import os
import pathlib
import pickle as pkl
import time

import lgdo.lh5 as lh5
import numpy as np
from legendmeta import LegendMetadata
from legendmeta.catalog import Props
from lgdo import Array, Table
from pygama.pargen.dplms_ge_dict import dplms_ge_dict

os.environ["LGDO_CACHE"] = "false"
os.environ["LGDO_BOUNDSCHECK"] = "false"
os.environ["DSPEED_CACHE"] = "false"
os.environ["DSPEED_BOUNDSCHECK"] = "false"
os.environ["PYGAMA_PARALLEL"] = "false"
os.environ["PYGAMA_FASTMATH"] = "false"

argparser = argparse.ArgumentParser()
argparser.add_argument("--fft_raw_filelist", help="fft_raw_filelist", type=str)
argparser.add_argument("--peak_file", help="tcm_filelist", type=str, required=True)
argparser.add_argument("--inplots", help="in_plot_path", type=str)

argparser.add_argument("--log", help="log_file", type=str)
argparser.add_argument("--database", help="database", type=str, required=True)
argparser.add_argument("--configs", help="configs", type=str, required=True)

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
logging.getLogger("lgdo").setLevel(logging.INFO)
logging.getLogger("h5py").setLevel(logging.INFO)
logging.getLogger("matplotlib").setLevel(logging.INFO)
logging.getLogger("dspeed.processing_chain").setLevel(logging.INFO)
logging.getLogger("legendmeta").setLevel(logging.INFO)

log = logging.getLogger(__name__)
sto = lh5.LH5Store()

configs = LegendMetadata(args.configs, lazy=True).on(args.timestamp, system=args.datatype)
dsp_config = configs["snakemake_rules"]["pars_dsp_dplms"]["inputs"]["proc_chain"][args.channel]

dplms_json = configs["snakemake_rules"]["pars_dsp_dplms"]["inputs"]["dplms_pars"][args.channel]
dplms_dict = Props.read_from(dplms_json)

db_dict = Props.read_from(args.database)

if dplms_dict["run_dplms"] is True:
    with open(args.fft_raw_filelist) as f:
        fft_files = sorted(f.read().splitlines())

    t0 = time.time()
    log.info("\nLoad fft data")
    energies = sto.read(f"{args.channel}/raw/daqenergy", fft_files)[0]
    idxs = np.where(energies.nda == 0)[0]
    raw_fft = sto.read(
        f"{args.channel}/raw", fft_files, n_rows=dplms_dict["n_baselines"], idx=idxs
    )[0]
    t1 = time.time()
    log.info(f"Time to load fft data {(t1-t0):.2f} s, total events {len(raw_fft)}")

    log.info("\nRunning event selection")
    peaks_kev = np.array(dplms_dict["peaks_kev"])
    kev_widths = [tuple(kev_width) for kev_width in dplms_dict["kev_widths"]]

    peaks_rounded = [int(peak) for peak in peaks_kev]
    peaks = sto.read(f"{args.channel}/raw", args.peak_file, field_mask=["peak"])[0]["peak"].nda
    ids = np.in1d(peaks, peaks_rounded)
    peaks = peaks[ids]
    idx_list = [np.where(peaks == peak)[0] for peak in peaks_rounded]

    raw_cal = sto.read(f"{args.channel}/raw", args.peak_file, idx=ids)[0]
    log.info(f"Time to run event selection {(time.time()-t1):.2f} s, total events {len(raw_cal)}")

    if isinstance(dsp_config, (str, list)):
        dsp_config = Props.read_from(dsp_config)

    if args.plot_path:
        out_dict, plot_dict = dplms_ge_dict(
            raw_fft,
            raw_cal,
            dsp_config,
            db_dict,
            dplms_dict,
            display=1,
        )
        if args.inplots:
            with open(args.inplots, "rb") as r:
                inplot_dict = pkl.load(r)
            inplot_dict.update({"dplms": plot_dict})

    else:
        out_dict = dplms_ge_dict(
            raw_fft,
            raw_cal,
            dsp_config,
            db_dict,
            dplms_dict,
        )

    coeffs = out_dict["dplms"].pop("coefficients")
    dplms_pars = Table(col_dict={"coefficients": Array(coeffs)})
    out_dict["dplms"][
        "coefficients"
    ] = f"loadlh5('{args.lh5_path}', '{args.channel}/dplms/coefficients')"

    log.info(f"DPLMS creation finished in {(time.time()-t0)/60} minutes")
else:
    out_dict = {}
    dplms_pars = Table(col_dict={"coefficients": Array([])})
    if args.inplots:
        with open(args.inplots, "rb") as r:
            inplot_dict = pkl.load(r)
    else:
        inplot_dict = {}

db_dict.update(out_dict)

pathlib.Path(os.path.dirname(args.lh5_path)).mkdir(parents=True, exist_ok=True)
sto.write(
    Table(col_dict={"dplms": dplms_pars}),
    name=args.channel,
    lh5_file=args.lh5_path,
    wo_mode="overwrite",
)

pathlib.Path(os.path.dirname(args.dsp_pars)).mkdir(parents=True, exist_ok=True)
Props.write_to(args.dsp_pars, dplms_dict)

if args.plot_path:
    pathlib.Path(os.path.dirname(args.plot_path)).mkdir(parents=True, exist_ok=True)
    with open(args.plot_path, "wb") as f:
        pkl.dump(inplot_dict, f, protocol=pkl.HIGHEST_PROTOCOL)
