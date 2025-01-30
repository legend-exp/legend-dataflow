import argparse
import logging
import pickle as pkl
import time
from pathlib import Path

import lgdo.lh5 as lh5
import numpy as np
from dbetto import TextDB
from dbetto.catalog import Props
from legendmeta import LegendMetadata
from lgdo import Array, Table
from pygama.pargen.dplms_ge_dict import dplms_ge_dict

from ..log import build_log

argparser = argparse.ArgumentParser()
argparser.add_argument("--fft_raw_filelist", help="fft_raw_filelist", type=str)
argparser.add_argument("--peak_file", help="tcm_filelist", type=str, required=True)
argparser.add_argument("--inplots", help="in_plot_path", type=str)
argparser.add_argument("--database", help="database", type=str, required=True)

argparser.add_argument("--log", help="log_file", type=str)
argparser.add_argument("--configs", help="configs", type=str, required=True)
argparser.add_argument("--metadata", help="metadata", type=str, required=True)

argparser.add_argument("--datatype", help="Datatype", type=str, required=True)
argparser.add_argument("--timestamp", help="Timestamp", type=str, required=True)
argparser.add_argument("--channel", help="Channel", type=str, required=True)

argparser.add_argument("--dsp_pars", help="dsp_pars", type=str, required=True)
argparser.add_argument("--lh5_path", help="lh5_path", type=str, required=True)
argparser.add_argument("--plot_path", help="plot_path", type=str)

args = argparser.parse_args()

configs = TextDB(args.configs, lazy=True).on(args.timestamp, system=args.datatype)
config_dict = configs["snakemake_rules"]["pars_dsp_dplms"]

log = build_log(config_dict, args.log)

log = logging.getLogger(__name__)
sto = lh5.LH5Store()

meta = LegendMetadata(path=args.metadata)
channel_dict = meta.channelmap(args.timestamp, system=args.datatype)
channel = f"ch{channel_dict[args.channel].daq.rawid:07}"

configs = LegendMetadata(args.configs, lazy=True).on(
    args.timestamp, system=args.datatype
)
dsp_config = config_dict["inputs"]["proc_chain"][args.channel]

dplms_json = config_dict["inputs"]["dplms_pars"][args.channel]
dplms_dict = Props.read_from(dplms_json)

db_dict = Props.read_from(args.database)

if dplms_dict["run_dplms"] is True:
    with Path(args.fft_raw_filelist).open() as f:
        fft_files = sorted(f.read().splitlines())

    t0 = time.time()
    log.info("\nLoad fft data")
    energies = sto.read(f"{channel}/raw/daqenergy", fft_files)[0]
    idxs = np.where(energies.nda == 0)[0]
    raw_fft = sto.read(
        f"{channel}/raw", fft_files, n_rows=dplms_dict["n_baselines"], idx=idxs
    )[0]
    t1 = time.time()
    log.info(f"Time to load fft data {(t1-t0):.2f} s, total events {len(raw_fft)}")

    log.info("\nRunning event selection")
    peaks_kev = np.array(dplms_dict["peaks_kev"])
    kev_widths = [tuple(kev_width) for kev_width in dplms_dict["kev_widths"]]

    peaks_rounded = [int(peak) for peak in peaks_kev]
    peaks = sto.read(f"{channel}/raw", args.peak_file, field_mask=["peak"])[0][
        "peak"
    ].nda
    ids = np.isin(peaks, peaks_rounded)
    peaks = peaks[ids]
    idx_list = [np.where(peaks == peak)[0] for peak in peaks_rounded]

    raw_cal = sto.read(f"{channel}/raw", args.peak_file, idx=ids)[0]
    log.info(
        f"Time to run event selection {(time.time()-t1):.2f} s, total events {len(raw_cal)}"
    )

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
            with Path(args.inplots).open("rb") as r:
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
    out_dict["dplms"]["coefficients"] = (
        f"loadlh5('{args.lh5_path}', '{channel}/dplms/coefficients')"
    )

    log.info(f"DPLMS creation finished in {(time.time()-t0)/60} minutes")
else:
    out_dict = {}
    dplms_pars = Table(col_dict={"coefficients": Array([])})
    if args.inplots:
        with Path(args.inplots).open("rb") as r:
            inplot_dict = pkl.load(r)
    else:
        inplot_dict = {}

db_dict.update(out_dict)

Path(args.lh5_path).parent.mkdir(parents=True, exist_ok=True)
sto.write(
    Table(col_dict={"dplms": dplms_pars}),
    name=channel,
    lh5_file=args.lh5_path,
    wo_mode="overwrite",
)

Path(args.dsp_pars).parent.mkdir(parents=True, exist_ok=True)
Props.write_to(args.dsp_pars, db_dict)

if args.plot_path:
    Path(args.plot_path).parent.mkdir(parents=True, exist_ok=True)
    with Path(args.plot_path).open("wb") as f:
        pkl.dump(inplot_dict, f, protocol=pkl.HIGHEST_PROTOCOL)
