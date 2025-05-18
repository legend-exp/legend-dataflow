from __future__ import annotations

import argparse
import pickle as pkl
import warnings
from pathlib import Path

import numpy as np
from dbetto import TextDB
from dbetto.catalog import Props
from legenddataflow import build_log
from legenddataflow.par.geds.pht.qc import build_qc

warnings.filterwarnings(action="ignore", category=RuntimeWarning)


def par_geds_pht_qc() -> None:
    argparser = argparse.ArgumentParser()
    argparser.add_argument("--cal-files", help="cal_files", nargs="*", type=str)
    argparser.add_argument("--fft-files", help="fft_files", nargs="*", type=str)

    argparser.add_argument(
        "--pulser-files", help="pulser_file", nargs="*", type=str, required=False
    )
    argparser.add_argument(
        "--overwrite-files", help="overwrite_files", nargs="*", type=str, required=False
    )

    argparser.add_argument("--configs", help="config", type=str, required=True)
    argparser.add_argument("--metadata", help="metadata path", type=str, required=True)
    argparser.add_argument("--log", help="log_file", type=str)

    argparser.add_argument("--datatype", help="Datatype", type=str, required=True)
    argparser.add_argument("--timestamp", help="Timestamp", type=str, required=True)
    argparser.add_argument("--channel", help="Channel", type=str, required=True)
    argparser.add_argument("--table-name", help="table name", type=str, required=True)

    argparser.add_argument(
        "--plot-path", help="plot_path", type=str, nargs="*", required=False
    )
    argparser.add_argument(
        "--save-path",
        help="save_path",
        type=str,
        nargs="*",
    )
    args = argparser.parse_args()

    configs = TextDB(args.configs, lazy=True).on(args.timestamp, system=args.datatype)
    config_dict = configs["snakemake_rules"]["pars_pht_qc"]

    build_log(config_dict, args.log)

    # get metadata dictionary
    channel_dict = config_dict["inputs"]["qc_config"][args.channel]
    kwarg_dict = Props.read_from(channel_dict)

    if args.overwrite_files:
        overwrite = Props.read_from(args.overwrite_files)
        if args.channel in overwrite:
            overwrite = overwrite[args.channel]["pars"]["operations"]
        else:
            overwrite = None
    else:
        overwrite = None

    # sort files in dictionary where keys are first timestamp from run
    if isinstance(args.cal_files, list):
        cal_files = []
        for file in args.cal_files:
            with Path(file).open() as f:
                cal_files += f.read().splitlines()
    else:
        with Path(args.cal_files).open() as f:
            cal_files = f.read().splitlines()

    cal_files = sorted(
        np.unique(cal_files)
    )  # need this as sometimes files get double counted as it somehow puts in the p%-* filelist and individual runs also

    if isinstance(args.fft_files, list):
        fft_files = []
        for file in args.fft_files:
            with Path(file).open() as f:
                fft_files += f.read().splitlines()
    else:
        with Path(args.fft_files).open() as f:
            fft_files = f.read().splitlines()

    fft_files = sorted(
        np.unique(fft_files)
    )  # need this as sometimes files get double counted as it somehow puts in the p%-* filelist and individual runs also

    hit_dict, plot_dict = build_qc(
        config=kwarg_dict,
        cal_files=cal_files,
        fft_files=fft_files,
        table_name=args.table_name,
        overwrite=overwrite,
        pulser_file=args.pulser_file,
        build_plots=bool(args.plot_path),
    )

    for file in args.save_path:
        Path(file).parent.mkdir(parents=True, exist_ok=True)
        Props.write_to(file, hit_dict)

    if args.plot_path:
        for file in args.plot_path:
            Path(file).parent.mkdir(parents=True, exist_ok=True)
            with Path(file).open("wb") as f:
                pkl.dump({"qc": plot_dict}, f, protocol=pkl.HIGHEST_PROTOCOL)
