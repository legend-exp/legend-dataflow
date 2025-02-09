"""
This script applies a simple calibration to the daqenergy for all channels,
it does this using a peak search, matching the peaks to the given ones
and deriving a simple scaling relation from adc to keV.
"""

import argparse
import logging
import pickle as pkl
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from dbetto.catalog import Props
from lgdo import lh5
from pygama.pargen.energy_cal import HPGeCalibration

from ....table_name import get_table_name

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

    argparser.add_argument("-d", "--debug", help="debug_mode", action="store_true")
    args = argparser.parse_args()

    logging.basicConfig(level=logging.DEBUG, filename=args.log, filemode="w")
    logging.getLogger("numba").setLevel(logging.INFO)
    logging.getLogger("parse").setLevel(logging.INFO)
    logging.getLogger("lgdo").setLevel(logging.INFO)
    logging.getLogger("matplotlib").setLevel(logging.INFO)
    log = logging.getLogger(__name__)

    channel = get_table_name(args.meta, args.timestamp, args.datatype, args.channel)

    # peaks to search for
    peaks_keV = np.array(
        [238, 583.191, 727.330, 860.564, 1592.53, 1620.50, 2103.53, 2614.50]
    )

    E_uncal = lh5.read(f"{channel}/raw/daqenergy", sorted(args.files))[0].view_as("np")
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
    detected_peaks_locs, detected_peaks_keV, roughpars = hpge_cal.hpge_find_E_peaks(
        E_uncal
    )

    log.info(f"{len(detected_peaks_locs)} peaks found:")
    log.info("\t   Energy   | Position  ")
    for i, (Li, Ei) in enumerate(zip(detected_peaks_locs, detected_peaks_keV)):
        log.info(f"\t{i}".ljust(4) + str(Ei).ljust(9) + f"| {Li:g}".ljust(5))  # noqa: G003

    # dictionary to pass to build hit
    out_dict = {
        "pars": {
            "operations": {
                "daqenergy_cal": {
                    "expression": "daqenergy*a",
                    "parameters": {"a": round(roughpars[0], 5)},
                }
            }
        }
    }

    # plot to check thagt the calibration is correct with zoom on 2.6 peak
    fig = plt.figure(figsize=(8, 10))
    ax = plt.subplot(211)
    ax.hist(E_uncal * roughpars[0], bins=np.arange(0, 3000, 1), histtype="step")
    ax.set_ylabel("counts")
    ax.set_yscale("log")
    ax2 = plt.subplot(212)
    ax2.hist(
        E_uncal * roughpars[0],
        bins=np.arange(2600, 2630, 1 * roughpars[0]),
        histtype="step",
    )
    ax2.set_xlabel("energy (keV)")
    ax2.set_ylabel("counts")
    plt.suptitle(args.channel)
    with Path(args.plot_file).open("wb") as w:
        pkl.dump(fig, w, protocol=pkl.HIGHEST_PROTOCOL)
    plt.close()

    Props.write_to_file(args.blind_curve, out_dict)
