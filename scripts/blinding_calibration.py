"""
This script applies a simple calibration to the daqenergy for all channels,
it does this using a peak search, matching the peaks to the given ones
and deriving a simple scaling relation from adc to keV.
"""

import argparse
import glob
import json
import logging
import os

import lgdo.lh5_store as lh5
import matplotlib.pyplot as plt
import numpy as np
from legendmeta import LegendMetadata
from matplotlib.backends.backend_pdf import PdfPages
from pygama.math.histogram import better_int_binning, get_hist
from pygama.pargen.energy_cal import hpge_find_E_peaks

argparser = argparse.ArgumentParser()
argparser.add_argument("input", help="input file path", type=str)
argparser.add_argument("output", help="output file", type=str)
argparser.add_argument("meta", help="meta path", type=str)
argparser.add_argument("out_plot", help="out plot path", type=str)
args = argparser.parse_args()


logging.basicConfig(level=logging.DEBUG)
logging.getLogger("numba").setLevel(logging.INFO)
logging.getLogger("parse").setLevel(logging.INFO)
logging.getLogger("lgdo").setLevel(logging.INFO)
logging.getLogger("matplotlib").setLevel(logging.INFO)

# get all files in specified directory
files = sorted(glob.glob(os.path.join(args.input, "*.lh5")))

# load in channel map
meta = LegendMetadata(args.meta)
chmap = meta.channelmap(os.path.basename(files[0]).split("-")[4])

# get list of all processable detectors
working_dets = [
    det
    for det, dic in chmap.map("system", unique=False)["geds"].map("name").items()
    if dic["analysis"]["processable"] is True
]

pars_dict = {}
# peaks to search for
peaks_keV = np.array([238, 583.191, 727.330, 860.564, 1592.53, 1620.50, 2103.53, 2614.50])

# get calobration for each channel
with PdfPages(args.out_plot) as pdf:
    for det in working_dets:
        E_uncal = lh5.load_nda(files, ["daqenergy"], f"ch{chmap[det].daq.rawid}/raw")["daqenergy"]

        guess_keV = 2620 / np.nanpercentile(E_uncal, 99)  # usual simple guess
        Euc_min = peaks_keV[0] / guess_keV * 0.6
        Euc_max = peaks_keV[-1] / guess_keV * 1.1
        dEuc = 1 / guess_keV

        # daqenergy is an int so use integer binning (dx used to be bugged as output so switched to nbins)
        Euc_min, Euc_max, nbins = better_int_binning(
            x_lo=Euc_min, x_hi=Euc_max, n_bins=(Euc_max - Euc_min) / dEuc
        )
        hist, bins, var = get_hist(E_uncal, range=(Euc_min, Euc_max), bins=nbins)

        # Run the rough peak search
        detected_peaks_locs, detected_peaks_keV, roughpars = hpge_find_E_peaks(
            hist, bins, var, peaks_keV, n_sigma=5, deg=0
        )
        # dictionary to pass to build hit
        out_dict = {
            "daqenergy_cal": {
                "expression": "daqenergy*a",
                "parameters": {"a": round(roughpars[0], 5)},
            }
        }
        pars_dict[f"ch{chmap[det].daq.rawid}"] = out_dict

        # plot to check thagt the calibration is correct with zoom on 2.6 peak
        fig = plt.figure()
        plt.hist(E_uncal * roughpars[0], bins=np.arange(0, 3000, 1), histtype="step")
        plt.xlabel("energy (keV)")
        plt.ylabel("counts")
        plt.yscale("log")
        plt.title(det)
        pdf.savefig()
        plt.close()


with open(args.output, "w") as w:
    json.dump(pars_dict, w, indent=4)
