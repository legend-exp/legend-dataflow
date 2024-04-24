import argparse
import json
import logging
import os
import pathlib

argparser = argparse.ArgumentParser()
argparser.add_argument("--log", help="log file", type=str)
argparser.add_argument("--output_file", help="output par file", type=str, required=True)
argparser.add_argument("--input_file", help="input par file", type=str, required=True)
argparser.add_argument("--svm_file", help="svm file", required=True)
args = argparser.parse_args()


if args.log is not None:
    pathlib.Path(os.path.dirname(args.log)).mkdir(parents=True, exist_ok=True)
    logging.basicConfig(level=logging.DEBUG, filename=args.log, filemode="w")
else:
    logging.basicConfig(level=logging.DEBUG)

logging.getLogger("parse").setLevel(logging.INFO)
logging.getLogger("lgdo").setLevel(logging.INFO)
logging.getLogger("h5py").setLevel(logging.INFO)

log = logging.getLogger(__name__)

with open(args.input_file) as r:
    par_data = json.load(r)

file = f"'$_/{os.path.basename(args.svm_file)}'"

par_data["svm"] = {"model_file": file}

pathlib.Path(os.path.dirname(args.output_file)).mkdir(parents=True, exist_ok=True)
with open(args.output_file, "w") as w:
    json.dump(par_data, w, indent=4)
