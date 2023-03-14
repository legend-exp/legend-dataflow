from datetime import datetime
from pygama.flow.file_db import FileDB
import argparse
import pathlib, os, json
from pathlib import Path


def check_log_files(log_path, output_file, gen_output):
    os.makedirs(os.path.dirname(output_file),exist_ok=True)
    with open(output_file, "w") as f:
        n_errors=0
        for file in Path(log_path).rglob("*.log"):
            with open(file) as r:
                text = r.read()
                if "ERROR" in text:
                    for line in text.splitlines():
                        if "ERROR" in line:
                            if n_errors ==0:
                                f.write(f"{gen_output} succesfully generated at {datetime.utcnow().strftime('%d/%m/%y %H:%M')} with errors \n")
                            f.write(f"{os.path.basename(file)} : {line}\n")
                else:
                    pass
            os.remove(file)
            text=None
        if n_errors ==0:
            f.write(f"{gen_output} succesfully generated at {datetime.utcnow().strftime('%d/%m/%y %H:%M')} with no errors \n")
    walk = list(os.walk(log_path))
    for path, _, _ in walk[::-1]:
        if len(os.listdir(path)) == 0:
            os.rmdir(path)

argparser = argparse.ArgumentParser()
argparser.add_argument("--log_path", help="tmp log path", type=str, required=True)
argparser.add_argument("--fileDBconfig", help="fileDBconfig", type=dict, required=False)
argparser.add_argument("--filelist", help="filelist", required=True, nargs='*',type=str)
argparser.add_argument("--gen_output", help="gen_output", type=str)
argparser.add_argument("--summary_log", help="summary_log", type=str)
argparser.add_argument("--fileDB", help="fileDB", type=str, required=False)
args = argparser.parse_args()

check_log_files(args.log_path, args.summary_log, args.gen_output)

# fdb = FileDB(args.fileDBconfig)
# fdb.scan_tables_columns()
# fdb.to_disk(args.fileDB)

pathlib.Path(args.gen_output).touch()