import argparse
import logging
import os
import pathlib

import numpy as np
from legendmeta import LegendMetadata
from legendmeta.catalog import Props
from pygama.skm.build_skm import build_skm
from util.FileKey import ProcessingFileKey


def group_files(fs_evt, fs_hit, fs_dsp, fs_tcm):
    """
    This function makes sure the files are ordered properly by matching the keys together
    returns list of tuples of (f_evt f_hit f_dsp f_tcm)
    """
    grouped_files = []
    for f_evt in sorted(
        fs_evt,
        key=lambda filename: ProcessingFileKey.get_filekey_from_pattern(
            os.path.basename(filename)
        ).get_unix_timestamp(),
    ):
        key = ProcessingFileKey.get_filekey_from_pattern(os.path.basename(f_evt))
        for f_dsp in fs_dsp:
            dsp_key = ProcessingFileKey.get_filekey_from_pattern(os.path.basename(f_dsp))
            if dsp_key.timestamp == key.timestamp:
                break
        for f_hit in fs_hit:
            hit_key = ProcessingFileKey.get_filekey_from_pattern(os.path.basename(f_hit))
            if hit_key.timestamp == key.timestamp:
                break
        for f_tcm in fs_tcm:
            tcm_key = ProcessingFileKey.get_filekey_from_pattern(os.path.basename(f_tcm))
            if tcm_key.timestamp == key.timestamp:
                break
        grouped_files.append((f_evt, f_hit, f_dsp, f_tcm))

    return grouped_files


argparser = argparse.ArgumentParser()
argparser.add_argument("--hit_files", help="hit files", nargs="*", type=str)
argparser.add_argument("--dsp_files", help="dsp files", nargs="*", type=str)
argparser.add_argument("--tcm_files", help="tcm files", nargs="*", type=str)
argparser.add_argument("--evt_files", help="evt files", nargs="*", type=str)

argparser.add_argument("--configs", help="configs", type=str, required=True)
argparser.add_argument("--datatype", help="Datatype", type=str, required=True)
argparser.add_argument("--timestamp", help="Timestamp", type=str, required=True)

argparser.add_argument("--metadata", help="metadata path", type=str, required=True)

argparser.add_argument("--log", help="log_file", type=str)

argparser.add_argument("--output", help="output file", type=str)
args = argparser.parse_args()

pathlib.Path(os.path.dirname(args.log)).mkdir(parents=True, exist_ok=True)
logging.basicConfig(level=logging.DEBUG, filename=args.log, filemode="w")
logging.getLogger("numba").setLevel(logging.INFO)
logging.getLogger("parse").setLevel(logging.INFO)
logging.getLogger("lgdo").setLevel(logging.INFO)
logging.getLogger("h5py._conv").setLevel(logging.INFO)

log = logging.getLogger(__name__)

# load in config
configs = LegendMetadata(path=args.configs)
skm_config_file = configs.on(args.timestamp, system=args.datatype)["snakemake_rules"]["tier_skm"][
    "inputs"
]["skm_config"]

skm_config = Props.read_from(skm_config_file)

if isinstance(args.hit_files, list) and args.hit_files[0].split(".")[-1] == "filelist":
    hit_files = args.hit_files[0]
    with open(hit_files) as f:
        hit_files = f.read().splitlines()
else:
    hit_files = args.hit_files

if isinstance(args.dsp_files, list) and args.dsp_files[0].split(".")[-1] == "filelist":
    dsp_files = args.dsp_files[0]
    with open(dsp_files) as f:
        dsp_files = f.read().splitlines()
else:
    dsp_files = args.dsp_files

if isinstance(args.evt_files, list) and args.evt_files[0].split(".")[-1] == "filelist":
    evt_files = args.evt_files[0]
    with open(evt_files) as f:
        evt_files = f.read().splitlines()
else:
    evt_files = args.evt_files

if isinstance(args.tcm_files, list) and args.tcm_files[0].split(".")[-1] == "filelist":
    tcm_files = args.tcm_files[0]
    with open(tcm_files) as f:
        tcm_files = f.read().splitlines()
else:
    tcm_files = args.tcm_files

if len(hit_files) != len(dsp_files) != len(tcm_files) != len(evt_files):
    err_msg = "Number of files for each tier must match"
    raise RuntimeError(err_msg)

log.debug(
    f"{len(hit_files)} hit files, {len(dsp_files)} dsp files, {len(tcm_files)} tcm files, {len(evt_files)} evt files"
)

input_files = group_files(evt_files, hit_files, dsp_files, tcm_files)

log.debug(f"{len(input_files)} file pairs")

pathlib.Path(os.path.dirname(args.output)).mkdir(parents=True, exist_ok=True)

rng = np.random.default_rng()
rand_num = f"{rng.integers(0,99999):05d}"
temp_output = f"{args.output}.{rand_num}"

for f_evt, f_hit, f_dsp, f_tcm in input_files:
    log_string = f"running files evt:{os.path.basename(f_evt)}, hit:{os.path.basename(f_hit)},"
    log_string += f"\ndsp:{os.path.basename(f_dsp)}, tcm: {os.path.basename(f_tcm)}"
    log.debug(log_string)
    build_skm(
        f_evt,
        f_hit,
        f_dsp,
        f_tcm,
        temp_output,
        skm_config,
        wo_mode="a",
        skm_group="skm",
        evt_group="evt",
        tcm_group="hardware_tcm_1",
        dsp_group="dsp",
        hit_group="hit",
        tcm_id_table_pattern="ch{}",
    )


os.rename(temp_output, args.output)
