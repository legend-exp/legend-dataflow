"""
This script applies a simple calibration to the daqenergy for all channels,
it does this using a peak search, matching the peaks to the given ones
and deriving a simple scaling relation from adc to keV.
"""

import argparse
import logging
import pickle as pkl

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from legendmeta import LegendMetadata
from legendmeta.catalog import Props
from lgdo import lh5
from pygama.math.histogram import better_int_binning, get_hist
from pygama.pargen.energy_cal import hpge_find_E_peaks

mpl.use("agg")

argparser = argparse.ArgumentParser()
argparser.add_argument("--files", help="files", nargs="*", type=str)
argparser.add_argument("--blind_curve", help="blind_curve", type=str)
argparser.add_argument("--plot_file", help="out plot path", type=str)
argparser.add_argument("--meta", help="meta", type=str)
argparser.add_argument("--timestamp", help="timestamp", type=str)
argparser.add_argument("--datatype", help="datatype", type=str)
argparser.add_argument("--channel", help="channel", type=str)
argparser.add_argument("--configs", help="configs", type=str)
argparser.add_argument("--log", help="log", type=str)
args = argparser.parse_args()

logging.basicConfig(level=logging.DEBUG, filename=args.log, filemode="w")
logging.getLogger("numba").setLevel(logging.INFO)
logging.getLogger("parse").setLevel(logging.INFO)
logging.getLogger("lgdo").setLevel(logging.INFO)
logging.getLogger("matplotlib").setLevel(logging.INFO)
log = logging.getLogger(__name__)

# load in channel map
meta = LegendMetadata(args.meta, lazy=True)
chmap = meta.channelmap(args.timestamp)

# if chmap.map("daq.rawid")[int(args.channel[2:])]["analysis"]["is_blinded"] is True:
pars_dict = {}
# peaks to search for
peaks_keV = np.array([238, 583.191, 727.330, 860.564, 1592.53, 1620.50, 2103.53, 2614.50])

E_uncal = lh5.read(f"{args.channel}/raw/daqenergy", sorted(args.files))[0].view_as("np")
E_uncal = E_uncal[E_uncal > 200]
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
ax2.hist(E_uncal * roughpars[0], bins=np.arange(2600, 2630, 1 * roughpars[0]), histtype="step")
ax2.set_xlabel("energy (keV)")
ax2.set_ylabel("counts")
plt.suptitle(args.channel)
with open(args.plot_file, "wb") as w:
    pkl.dump(fig, w, protocol=pkl.HIGHEST_PROTOCOL)
plt.close()

# else:
#     out_dict = {
#         "pars": {
#             "operations": {
#                 "daqenergy_cal": {
#                     "expression": "daqenergy*a",
#                     "parameters": {"a": np.nan},
#                 }
#             }
#         }
#     }
#     fig = plt.figure(figsize=(8, 10))
#     plt.suptitle(f"{args.channel}-blind_off")
#     with open(args.plot_file, "wb") as w:
#         pkl.dump(fig, w, protocol=pkl.HIGHEST_PROTOCOL)
#     plt.close()
Props.write_to_file(args.blind_curve, out_dict)
