import argparse

from legendmeta.catalog import Props
from pygama.flow.file_db import FileDB

argparser = argparse.ArgumentParser()
argparser.add_argument("--config", required=True)
argparser.add_argument("--scan-path", required=True)
argparser.add_argument("--output", required=True)
args = argparser.parse_args()

config = Props.read_from(args.config)

fdb = FileDB(config, scan=False)
fdb.scan_files([args.scan_path])
fdb.scan_tables_columns(dir_files_conform=True)
fdb.to_disk(args.output, wo_mode="of")
