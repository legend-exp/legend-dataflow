
import argparse
import logging
import os
import pickle as pkl
import lgdo.lh5 as lh5
import lgdo.types as lgd
import numpy as np
from legendmeta.catalog import Props

os.environ["LGDO_CACHE"] = "false"
os.environ["LGDO_BOUNDSCHECK"] = "false"

argparser = argparse.ArgumentParser()
argparser.add_argument("--log", help="log file", type=str)
argparser.add_argument("--output", help="output xtc matrix", type=str, required=True)
argparser.add_argument("--xtc_matrix", help="input xtc matrix", type=str, required=True)
argparser.add_argument("--inpars", help="input cal pars", required=True, nargs="*")
args = argparser.parse_args()

logging.basicConfig(level=logging.DEBUG, filename=args.log, filemode="w")
logging.getLogger("parse").setLevel(logging.INFO)
logging.getLogger("lgdo").setLevel(logging.INFO)
logging.getLogger("h5py").setLevel(logging.INFO)

sto = lh5.LH5Store()
log = logging.getLogger(__name__)

# Load files
xtc, _ = sto.read("xtc", args.xtc_matrix)
log.debug("loaded matrix")

par_pht = Props.read_from(
    args.inpars
)

gains = {}
for chan in xtc.rawid_index.nda:
    try:
        gains[chan] = par_pht[f"ch{chan}"]["pars"]["operations"]["cuspEmax_ctc_runcal"]["parameters"]["b"]
    except:
        gains[chan] =np.nan
    
out_matrix_positive = np.full_like(xtc.xtalk_matrix_positive.nda, np.nan)
out_matrix_negative = np.full_like(xtc.xtalk_matrix_negative.nda, np.nan)
for idx, _ in np.ndenumerate(out_matrix_positive):
    i, j = idx
    gain_corr = gains[xtc.rawid_index.nda[j]]/gains[xtc.rawid_index.nda[i]]
    out_matrix_positive[i][j] = xtc.xtalk_matrix_positive.nda[i][j] * gain_corr
    out_matrix_negative[i][j] = xtc.xtalk_matrix_positive.nda[i][j] * gain_corr

out = lgd.Table(col_dict={"rawid_index": xtc.rawid_index, 
                        "xtalk_matrix_positive":lgd.Array(nda=out_matrix_positive),
                        "xtalk_matrix_negative": lgd.Array(nda=out_matrix_negative)})

sto.write(out, "xtc", lh5_file=args.output)