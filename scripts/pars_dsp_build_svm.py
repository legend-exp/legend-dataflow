import argparse
import logging
import os
import pickle as pkl

os.environ["LGDO_CACHE"] = "false"
os.environ["LGDO_BOUNDSCHECK"] = "false"
os.environ["DSPEED_CACHE"] = "false"
os.environ["DSPEED_BOUNDSCHECK"] = "false"

import lgdo.lh5 as lh5
from legendmeta.catalog import Props
from sklearn.svm import SVC

argparser = argparse.ArgumentParser()
argparser.add_argument("--log", help="log file", type=str)
argparser.add_argument("--output_file", help="output SVM file", type=str, required=True)
argparser.add_argument("--train_data", help="input data file", type=str, required=True)
argparser.add_argument("--train_hyperpars", help="input hyperparameter file", required=True)
args = argparser.parse_args()

logging.basicConfig(level=logging.DEBUG, filename=args.log, filemode="w")
logging.getLogger("parse").setLevel(logging.INFO)
logging.getLogger("lgdo").setLevel(logging.INFO)
logging.getLogger("h5py").setLevel(logging.INFO)

sto = lh5.LH5Store()
log = logging.getLogger(__name__)

# Load files
tb, _ = sto.read("ml_train/dsp", args.train_data)
log.debug("loaded data")

hyperpars = Props.read_from(args.train_hyperpars)

# Define training inputs
dwts_norm = tb["dwt_norm"].nda
labels = tb["dc_label"].nda


log.debug("training model")
# Initialize and train SVM
svm = SVC(
    random_state=int(hyperpars["random_state"]),
    kernel=hyperpars["kernel"],
    decision_function_shape=hyperpars["decision_function_shape"],
    class_weight=hyperpars["class_weight"],
    C=float(hyperpars["C"]),
    gamma=float(hyperpars["gamma"]),
)

svm.fit(dwts_norm, labels)

log.debug("trained model")

# Save trained model with pickle
with open(args.output_file, "wb") as svm_file:
    pkl.dump(svm, svm_file, protocol=pkl.HIGHEST_PROTOCOL)
