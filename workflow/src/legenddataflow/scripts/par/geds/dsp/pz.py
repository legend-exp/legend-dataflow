import argparse
import copy
import pickle as pkl
from pathlib import Path

import lgdo.lh5 as lh5
import numpy as np
from dbetto import TextDB
from dbetto.catalog import Props
from pygama.pargen.data_cleaning import get_cut_indexes
from pygama.pargen.dsp_optimize import run_one_dsp
from pygama.pargen.pz_correct import PZCorrect

from .....convert_np import convert_dict_np_to_float
from .....log import build_log
from ....pulser_removal import get_pulser_mask


def par_geds_dsp_pz() -> None:
    argparser = argparse.ArgumentParser()
    argparser.add_argument("--configs", help="configs path", type=str, required=True)
    argparser.add_argument("--log", help="log file", type=str)
    argparser.add_argument(
        "-p", "--no-pulse", help="no pulser present", action="store_true"
    )

    argparser.add_argument("--datatype", help="Datatype", type=str, required=True)
    argparser.add_argument("--timestamp", help="Timestamp", type=str, required=True)
    argparser.add_argument("--channel", help="Channel", type=str, required=True)
    argparser.add_argument(
        "--raw-table-name", help="raw table name", type=str, required=True
    )

    argparser.add_argument("--plot-path", help="plot path", type=str, required=False)
    argparser.add_argument("--output-file", help="output file", type=str, required=True)

    argparser.add_argument(
        "--pulser-file", help="pulser file", type=str, required=False
    )

    argparser.add_argument("--raw-files", help="input files", nargs="*", type=str)
    argparser.add_argument("--pz-files", help="input files", nargs="*", type=str)
    args = argparser.parse_args()

    configs = TextDB(args.configs, lazy=True).on(args.timestamp, system=args.datatype)
    config_dict = configs["snakemake_rules"]["pars_dsp_tau"]

    log = build_log(config_dict, args.log)

    channel_dict = config_dict["inputs"]["processing_chain"][args.channel]
    kwarg_dict = config_dict["inputs"]["tau_config"][args.channel]

    kwarg_dict = Props.read_from(kwarg_dict)

    if kwarg_dict["run_tau"] is True:
        dsp_config = Props.read_from(channel_dict)
        kwarg_dict.pop("run_tau")
        if args.pz_files is not None and len(args.pz_files) > 0:
            if (
                isinstance(args.pz_files, list)
                and args.pz_files[0].split(".")[-1] == "filelist"
            ):
                input_file = args.pz_files[0]
                with Path(input_file).open() as f:
                    input_file = f.read().splitlines()
            else:
                input_file = args.pz_files
        if len(input_file) == 0:
            if (
                isinstance(args.raw_files, list)
                and args.raw_files[0].split(".")[-1] == "filelist"
            ):
                input_file = args.raw_files[0]
                with Path(input_file).open() as f:
                    input_file = f.read().splitlines()
            else:
                input_file = args.raw_files

        log.debug(f"Reading Data for {args.raw_table_name} from:")
        log.debug(input_file)

        data = lh5.read(
            args.raw_table_name,
            input_file,
            field_mask=["daqenergy", "timestamp", "t_sat_lo"],
        ).view_as("pd")
        threshold = kwarg_dict.pop("threshold")

        if args.no_pulse is False and (
            args.pz_files is None or len(args.pz_files) == 0
        ):
            mask = get_pulser_mask(args.pulser_file)
        else:
            mask = np.full(len(data), False)

        discharges = data["t_sat_lo"] > 0
        discharge_timestamps = np.where(data["timestamp"][discharges])[0]
        is_recovering = np.full(len(data), False, dtype=bool)
        for tstamp in discharge_timestamps:
            is_recovering = is_recovering | np.where(
                (
                    ((data["timestamp"] - tstamp) < 0.01)
                    & ((data["timestamp"] - tstamp) > 0)
                ),
                True,
                False,
            )
        cuts = np.where(
            (data.daqenergy.to_numpy() > threshold) & (~mask) & (~is_recovering)
        )[0]
        log.debug(f"{len(cuts)} events passed threshold and pulser cuts")
        log.debug(cuts)
        tb_data = lh5.read(
            args.raw_table_name,
            input_file,
            idx=cuts,
            n_rows=kwarg_dict["n_events"] * 2,
        )

        dsp_config_optimise_removed = copy.deepcopy(dsp_config)
        if "tau1" in dsp_config["outputs"]:
            dsp_config_optimise_removed["outputs"].remove("tau1")
        if "tau2" in dsp_config["outputs"]:
            dsp_config_optimise_removed["outputs"].remove("tau2")
        if "frac" in dsp_config["outputs"]:
            dsp_config_optimise_removed["outputs"].remove("frac")

        tb_out = run_one_dsp(tb_data, dsp_config_optimise_removed)
        log.debug("Processed Data")
        cut_parameters = kwarg_dict.get("cut_parameters", None)
        if cut_parameters is not None:
            idxs = get_cut_indexes(tb_out, cut_parameters=cut_parameters)
            log.debug("Applied cuts")
            log.debug(f"{len(idxs)} events passed cuts")
            tb_data = lh5.read(
                args.raw_table_name,
                input_file,
                idx=cuts[: 2 * kwarg_dict["n_events"]][idxs],
                n_rows=kwarg_dict.pop("n_events"),
            )

        tau = PZCorrect(
            dsp_config,
            kwarg_dict["wf_field"],
            debug_mode=kwarg_dict.get("debug_mode", False),
        )
        log.debug("Calculating pz constant")
        if kwarg_dict["mode"] == "single":
            tau.get_single_decay_constant(
                tb_data, kwarg_dict.get("slope_param", "tail_slope")
            )
            log.debug(
                f"Found tau: {tau.output_dict['pz']['tau']}+- {tau.output_dict['pz']['tau_err']}"
            )
        elif kwarg_dict["mode"] == "double":
            tau.get_dpz_decay_constants(
                tb_data,
                kwarg_dict.get("percent_tau1_fit", 0.1),
                kwarg_dict.get("percent_tau2_fit", 0.2),
                kwarg_dict.get("offset_from_wf_max", 10),
                kwarg_dict.get("superpulse_bl_idx", 25),
                kwarg_dict.get("superpulse_window_width", 13),
            )
            log.debug("found dpz constants : ")
            for entry in ["tau1", "tau2", "frac"]:
                log.debug(
                    f"{entry}:{tau.output_dict['pz'][entry]}+- {tau.output_dict['pz'][f'{entry}_err']}"
                )
        else:
            msg = f"Unknown mode: {kwarg_dict['mode']}, must be either single or double"
            raise ValueError(msg)
        tau.dsp_config = dsp_config_optimise_removed

        if args.plot_path:
            Path(args.plot_path).parent.mkdir(parents=True, exist_ok=True)

            plot_dict = tau.plot_waveforms_after_correction(
                tb_data, "wf_pz", norm_param=kwarg_dict.get("norm_param", "pz_mean")
            )

            zoomed = tau.plot_waveforms_after_correction(
                tb_data,
                "wf_pz",
                norm_param=kwarg_dict.get("norm_param", "pz_mean"),
                xlim=[400, 1020],
                ylim=[0.8, 1.1],
            )

            plot_dict.update({"waveforms_zoomed": zoomed["waveforms"]})

            plot_dict.update(
                tau.plot_slopes(
                    tb_data, kwarg_dict.get("final_slope_param", "pz_slope")
                )
            )
            plot_dict.update(
                tau.plot_slopes(
                    tb_data, kwarg_dict.get("final_slope_param", "pz_slope"), True
                )
            )

            with Path(args.plot_path).open("wb") as f:
                pkl.dump({"tau": plot_dict}, f, protocol=pkl.HIGHEST_PROTOCOL)
        out_dict = convert_dict_np_to_float(tau.output_dict)
    else:
        out_dict = {}

    Path(args.output_file).parent.mkdir(parents=True, exist_ok=True)
    Props.write_to(args.output_file, out_dict)
