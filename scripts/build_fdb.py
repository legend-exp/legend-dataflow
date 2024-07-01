import argparse

from legendmeta.catalog import Props
from pygama.flow.file_db import FileDB

argparser = argparse.ArgumentParser()
argparser.add_argument("--config", help="config", type=str, required=True)
argparser.add_argument("--file_path", help="files_path", type=str, required=True)
argparser.add_argument("--output_file", help="output_file", type=str, required=True)
args = argparser.parse_args()

config = Props.read_from(args.config)

fdb = FileDB(config, scan=False)
fdb.scan_files([args.file_path])
fdb.scan_tables_columns(dir_files_conform=True)
fdb.to_disk(args.output_file, wo_mode="of")
