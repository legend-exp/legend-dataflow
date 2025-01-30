import argparse
import pickle as pkl
from pathlib import Path

from dbetto import TextDB
from dbetto.catalog import Props
from lgdo import lh5
from sklearn.svm import SVC

from ..log import build_log

argparser = argparse.ArgumentParser()
argparser.add_argument("--log", help="log file", type=str)
argparser.add_argument("--configs", help="config file", type=str)

argparser.add_argument("--datatype", help="Datatype", type=str, required=True)
argparser.add_argument("--timestamp", help="Timestamp", type=str, required=True)

argparser.add_argument("--output_file", help="output SVM file", type=str, required=True)
argparser.add_argument("--train_data", help="input data file", type=str, required=True)
argparser.add_argument("--train_hyperpars", help="input hyperparameter file", required=True)
args = argparser.parse_args()

configs = TextDB(args.configs, lazy=True).on(args.timestamp, system=args.datatype)
config_dict = configs["snakemake_rules"]["pars_dsp_build_svm"]

log = build_log(config_dict, args.log)

# Load files
tb = lh5.read("ml_train/dsp", args.train_data)
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
with Path(args.output_file).open("wb") as svm_file:
    pkl.dump(svm, svm_file, protocol=pkl.HIGHEST_PROTOCOL)
