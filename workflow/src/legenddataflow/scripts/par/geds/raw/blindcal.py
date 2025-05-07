"""
This script applies a simple calibration to the daqenergy for all channels,
it does this using a peak search, matching the peaks to the given ones
and deriving a simple scaling relation from adc to keV.
"""

from __future__ import annotations

import argparse
import pickle as pkl
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from dbetto import TextDB
from dbetto.catalog import Props
from legenddataflowscripts.utils import build_log
from lgdo import lh5
from pygama.pargen.energy_cal import HPGeCalibration

mpl.use("agg")


def par_geds_raw_blindcal() -> None:
    argparser = argparse.ArgumentParser()
    argparser.add_argument("--files", help="files", nargs="*", type=str)

    argparser.add_argument("--blind-curve", help="blind_curve", type=str)
    argparser.add_argument("--plot-file", help="out plot path", type=str)

    argparser.add_argument("--meta", help="meta", type=str)
    argparser.add_argument("--configs", help="configs", type=str)
    argparser.add_argument("--log", help="log", type=str)

    argparser.add_argument("--timestamp", help="timestamp", type=str)
    argparser.add_argument("--datatype", help="datatype", type=str)
    argparser.add_argument("--channel", help="channel", type=str)
    argparser.add_argument(
        "--raw-table-name", help="raw table name", type=str, required=True
    )

    argparser.add_argument("-d", "--debug", help="debug_mode", action="store_true")
    args = argparser.parse_args()

    configs = TextDB(args.configs, lazy=True).on(args.timestamp, system=args.datatype)
    config_dict = configs["snakemake_rules"]["tier_raw_blind_check"]

    log = build_log(config_dict, args.log)

    # peaks to search for
    peaks_keV = np.array(
        [238, 583.191, 727.330, 860.564, 1592.53, 1620.50, 2103.53, 2614.50]
    )

    if isinstance(args.files, list) and args.files[0].split(".")[-1] == "filelist":
        input_file = args.files[0]
        with Path(input_file).open() as f:
            input_file = f.read().splitlines()
    else:
        input_file = args.files

    E_uncal = lh5.read_as(f"{args.raw_table_name}/daqenergy", input_file, library="np")
    E_uncal = E_uncal[E_uncal > 200]
    guess_keV = 2620 / np.nanpercentile(E_uncal, 99)  # usual simple guess

    # daqenergy is an int so use integer binning (dx used to be bugged as output so switched to nbins)

    hpge_cal = HPGeCalibration(
        "daqenergy",
        peaks_keV,
        guess_keV,
        0,
        uncal_is_int=True,
        debug_mode=args.debug,
    )

    # Run the rough peak search
    hpge_cal.hpge_find_energy_peaks(E_uncal, etol_kev=5)
    detected_peaks_locs = hpge_cal.peak_locs
    detected_peaks_keV = hpge_cal.peaks_kev
    roughpars = hpge_cal.pars

    msg = f"{len(detected_peaks_locs)} peaks found:"
    log.info(msg)
    log.info("\t   Energy   | Position  ")
    for i, (Li, Ei) in enumerate(
        zip(detected_peaks_locs, detected_peaks_keV, strict=False)
    ):
        log.info(f"\t{i}".ljust(4) + str(Ei).ljust(9) + f"| {Li:g}".ljust(5))  # noqa: G003

    # dictionary to pass to build hit
    out_dict = {
        "pars": {
            "operations": {
                "daqenergy_cal": {
                    "expression": "daqenergy*a",
                    "parameters": {"a": round(float(roughpars[1]), 5)},
                }
            }
        }
    }

    # plot to check thagt the calibration is correct with zoom on 2.6 peak
    fig = plt.figure(figsize=(8, 10))
    ax = plt.subplot(211)
    ax.hist(E_uncal * roughpars[1], bins=np.arange(0, 3000, 1), histtype="step")
    ax.set_ylabel("counts")
    ax.set_yscale("log")
    ax2 = plt.subplot(212)
    ax2.hist(
        E_uncal * roughpars[1],
        bins=np.arange(2600, 2630, 1 * roughpars[1]),
        histtype="step",
    )
    ax2.set_xlabel("energy (keV)")
    ax2.set_ylabel("counts")
    plt.suptitle(args.channel)
    with Path(args.plot_file).open("wb") as w:
        pkl.dump(fig, w, protocol=pkl.HIGHEST_PROTOCOL)
    plt.close()

    Props.write_to(args.blind_curve, out_dict)
