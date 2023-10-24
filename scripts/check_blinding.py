import argparse
import logging
import os
import pathlib
from legendmeta.catalog import Props
import numpy as np
import numexpr as ne
import lgdo.lh5_store as lh5
from pygama.pargen.energy_cal import get_i_local_maxima
from pygama.math.histogram import get_hist
import matplotlib.pyplot as plt
import matplotlib as mpl 
mpl.use("Agg")

argparser = argparse.ArgumentParser()
argparser.add_argument("files", help="input file", type=str)
argparser.add_argument("output", help="output file", type=str)
argparser.add_argument("blind_curve", help="blinding curves file", type=str)
argparser.add_argument("--datatype", help="Datatype", type=str, required=True)
argparser.add_argument("--timestamp", help="Timestamp", type=str, required=True)
argparser.add_argument("--configs", help="config file", type=str)
argparser.add_argument("--channel", help="channel", type=str)
argparser.add_argument("--log", help="log file", type=str)
args = argparser.parse_args()

os.makedirs(os.path.dirname(args.log), exist_ok=True)
logging.basicConfig(level=logging.INFO, filename=args.log, filemode="w")

blind_curve = Props.read_from(args.blind_curve)[args.channel]

daqenergy = lh5.load_nda(sorted(args.files), ["daqenergy"],f"{channel}/raw")["daqenergy"]

# calibrate daq energy using pre existing curve
daqenergy_cal = ne.evaluate(
                blind_curve["daqenergy_cal"]["expression"],
                local_dict=dict(daqenergy=daqenergy, **blind_curve["daqenergy_cal"]["parameters"])
            )

# bin with 1 keV bins and get maxs
hist, bins, var = get_hist(daqenergy_cal, np.arange(0,3000,1))
maxs = get_i_local_maxima(hist, delta =5)
plt.figure()
plt.step((bins[1:]+bins[:-1])/2, hist, where="mid")
plt.close()


# check for peaks within +- 5keV of  2614 and 583 to ensure blinding still valid and if so create file else raise error
if np.any(np.abs(maxs-2614)<5) and np.any(np.abs(maxs-583)<5):
    pathlib.Path(os.path.dirname(args.output)).mkdir(parents=True, exist_ok=True)
    with open(args.output,"w") as f:
        json.dump({}, f)
else:
    raise RuntimeError("peaks not found in daqenergy")