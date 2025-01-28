import argparse
from pathlib import Path

from dbetto.catalog import Props

argparser = argparse.ArgumentParser()
argparser.add_argument("--log", help="log file", type=str)
argparser.add_argument("--output_file", help="output par file", type=str, required=True)
argparser.add_argument("--input_file", help="input par file", type=str, required=True)
argparser.add_argument("--svm_file", help="svm file", required=True)
args = argparser.parse_args()

par_data = Props.read_from(args.input_file)

file = f"'$_/{Path(args.svm_file).name}'"

par_data["svm"] = {"model_file": file}

Path(args.output_file).parent.mkdir(parents=True, exist_ok=True)
Props.write_to(args.output_file, par_data)
