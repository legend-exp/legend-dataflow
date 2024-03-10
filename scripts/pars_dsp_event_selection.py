import argparse
import json
import logging
import os
import pathlib
import pickle as pkl
import time
import warnings

os.environ["LGDO_CACHE"] = "false"
os.environ["LGDO_BOUNDSCHECK"] = "false"
os.environ["DSPEED_CACHE"] = "false"
os.environ["DSPEED_BOUNDSCHECK"] = "false"

import lgdo.lh5 as lh5
import numpy as np
import pygama.pargen.energy_optimisation as om
from legendmeta import LegendMetadata
from legendmeta.catalog import Props
from pygama.pargen.utils import get_tcm_pulser_ids

warnings.filterwarnings(action="ignore", category=RuntimeWarning)

argparser = argparse.ArgumentParser()
argparser.add_argument("--raw_filelist", help="raw_filelist", type=str)
argparser.add_argument("--tcm_filelist", help="tcm_filelist", type=str, required=True)
argparser.add_argument("--decay_const", help="decay_const", type=str, required=True)
argparser.add_argument("--configs", help="configs", type=str, required=True)

argparser.add_argument("--log", help="log_file", type=str)

argparser.add_argument("--datatype", help="Datatype", type=str, required=True)
argparser.add_argument("--timestamp", help="Timestamp", type=str, required=True)
argparser.add_argument("--channel", help="Channel", type=str, required=True)

argparser.add_argument("--peak_file", help="peak_file", type=str, required=True)
args = argparser.parse_args()

logging.basicConfig(level=logging.DEBUG, filename=args.log, filemode="w")
logging.getLogger("numba").setLevel(logging.INFO)
logging.getLogger("parse").setLevel(logging.INFO)
logging.getLogger("lgdo").setLevel(logging.INFO)
logging.getLogger("h5py").setLevel(logging.INFO)
logging.getLogger("matplotlib").setLevel(logging.INFO)
logging.getLogger("dspeed.processing_chain").setLevel(logging.INFO)


log = logging.getLogger(__name__)

t0 = time.time()

conf = LegendMetadata(path=args.configs)
configs = conf.on(args.timestamp, system=args.datatype)
dsp_config = configs["snakemake_rules"]["pars_dsp_peak_selection"]["inputs"]["processing_chain"][
    args.channel
]
peak_json = configs["snakemake_rules"]["pars_dsp_peak_selection"]["inputs"]["peak_config"][args.channel]

peak_dict = Props.read_from(opt_json)
db_dict = Props.read_from(args.decay_const)

if opt_dict.pop("run_selection") is True:
    with open(args.raw_filelist) as f:
        files = f.read().splitlines()

    raw_files = sorted(files)

    # get pulser mask from tcm files
    with open(args.tcm_filelist) as f:
        tcm_files = f.read().splitlines()
    tcm_files = sorted(np.unique(tcm_files))
    ids, mask = get_tcm_pulser_ids(
        tcm_files, args.channel, peak_dict["pulser_multiplicity_threshold"]
    )

    sto = lh5.LH5Store()
    idx_events, idx_list = om.event_selection(
        raw_files,
        f"{args.channel}/raw",
        dsp_config,
        db_dict,
        peaks_keV,
        np.arange(0, len(peaks_keV), 1).tolist(),
        kev_widths,
        pulser_mask=mask,
        cut_parameters=peak_dict["cut_parameters"],
        n_events=peak_dict["n_events"],
        threshold=peak_dict["threshold"],
        wf_field=peak_dict["wf_field"],
    )

    tb_data = sto.read(
        f"{args.channel}/raw",
        raw_files,
        idx=idx_events,
        n_rows=opt_dict["n_events"],
    )[0]

    pathlib.Path(os.path.dirname(args.peak_file)).mkdir(parents=True, exist_ok=True)
    sto.write(
        tb_data,
        name="raw",
        lh5_file=args.peak_file,
        wo_mode="overwrite",
    )
else:
    pathlib.Path(os.path.dirname(args.peak_file)).mkdir(parents=True, exist_ok=True)
    pathlib.Path(args.peak_file).touch()