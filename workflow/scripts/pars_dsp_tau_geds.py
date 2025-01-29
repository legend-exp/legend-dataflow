import argparse
import pickle as pkl
from pathlib import Path

import lgdo.lh5 as lh5
import numpy as np
from dbetto import TextDB
from dbetto.catalog import Props
from legenddataflow.log import build_log
from legendmeta import LegendMetadata
from pygama.pargen.data_cleaning import get_cut_indexes, get_tcm_pulser_ids
from pygama.pargen.dsp_optimize import run_one_dsp
from pygama.pargen.extract_tau import ExtractTau

argparser = argparse.ArgumentParser()
argparser.add_argument("--configs", help="configs path", type=str, required=True)
argparser.add_argument("--metadata", help="metadata", type=str, required=True)
argparser.add_argument("--log", help="log file", type=str)

argparser.add_argument("--datatype", help="Datatype", type=str, required=True)
argparser.add_argument("--timestamp", help="Timestamp", type=str, required=True)
argparser.add_argument("--channel", help="Channel", type=str, required=True)

argparser.add_argument("--plot_path", help="plot path", type=str, required=False)
argparser.add_argument("--output_file", help="output file", type=str, required=True)

argparser.add_argument("--pulser_file", help="pulser file", type=str, required=False)

argparser.add_argument("--raw_files", help="input files", nargs="*", type=str)
argparser.add_argument("--tcm_files", help="tcm_files", nargs="*", type=str, required=False)
args = argparser.parse_args()

sto = lh5.LH5Store()

configs = TextDB(args.configs, lazy=True).on(args.timestamp, system=args.datatype)
config_dict = configs["snakemake_rules"]["pars_dsp_nopt"]

log = build_log(config_dict, args.log)

meta = LegendMetadata(path=args.metadata)
channel_dict = meta.channelmap(args.timestamp, system=args.datatype)
channel = f"ch{channel_dict[args.channel].daq.rawid:07}"

channel_dict = config_dict["inputs"]["processing_chain"][args.channel]
kwarg_dict = config_dict["inputs"]["tau_config"][args.channel]

kwarg_dict = Props.read_from(kwarg_dict)

if kwarg_dict["run_tau"] is True:
    dsp_config = Props.read_from(channel_dict)
    kwarg_dict.pop("run_tau")
    if isinstance(args.raw_files, list) and args.raw_files[0].split(".")[-1] == "filelist":
        input_file = args.raw_files[0]
        with Path(input_file).open() as f:
            input_file = f.read().splitlines()
    else:
        input_file = args.raw_files

    if args.pulser_file:
        pulser_dict = Props.read_from(args.pulser_file)
        mask = np.array(pulser_dict["mask"])

    elif args.tcm_filelist:
        # get pulser mask from tcm files
        with Path(args.tcm_filelist).open() as f:
            tcm_files = f.read().splitlines()
        tcm_files = sorted(np.unique(tcm_files))
        ids, mask = get_tcm_pulser_ids(
            tcm_files, channel, kwarg_dict["pulser_multiplicity_threshold"]
        )
    else:
        msg = "No pulser file or tcm filelist provided"
        raise ValueError(msg)

    data = sto.read(
        f"{channel}/raw", input_file, field_mask=["daqenergy", "timestamp", "t_sat_lo"]
    )[0].view_as("pd")
    threshold = kwarg_dict.pop("threshold")

    discharges = data["t_sat_lo"] > 0
    discharge_timestamps = np.where(data["timestamp"][discharges])[0]
    is_recovering = np.full(len(data), False, dtype=bool)
    for tstamp in discharge_timestamps:
        is_recovering = is_recovering | np.where(
            (((data["timestamp"] - tstamp) < 0.01) & ((data["timestamp"] - tstamp) > 0)),
            True,
            False,
        )
    cuts = np.where((data.daqenergy.to_numpy() > threshold) & (~mask) & (~is_recovering))[0]

    tb_data = sto.read(
        f"{channel}/raw",
        input_file,
        idx=cuts,
        n_rows=kwarg_dict.pop("n_events"),
    )[0]

    tb_out = run_one_dsp(tb_data, dsp_config)
    log.debug("Processed Data")
    cut_parameters = kwarg_dict.get("cut_parameters", None)
    if cut_parameters is not None:
        idxs = get_cut_indexes(tb_out, cut_parameters=cut_parameters)
        log.debug("Applied cuts")
        log.debug(f"{len(idxs)} events passed cuts")
    else:
        idxs = np.full(len(tb_out), True, dtype=bool)

    tau = ExtractTau(dsp_config, kwarg_dict["wf_field"])
    slopes = tb_out["tail_slope"].nda
    log.debug("Calculating pz constant")

    tau.get_decay_constant(slopes[idxs], tb_data[kwarg_dict["wf_field"]])

    if args.plot_path:
        Path(args.plot_path).parent.mkdir(parents=True, exist_ok=True)

        plot_dict = tau.plot_waveforms_after_correction(
            tb_data, "wf_pz", norm_param=kwarg_dict.get("norm_param", "pz_mean")
        )
        plot_dict.update(tau.plot_slopes(slopes[idxs]))

        with Path(args.plot_path).open("wb") as f:
            pkl.dump({"tau": plot_dict}, f, protocol=pkl.HIGHEST_PROTOCOL)
else:
    out_dict = {}

Path(args.output_file).parent.mkdir(parents=True, exist_ok=True)
Props.write_to(args.output_file, tau.output_dict)
