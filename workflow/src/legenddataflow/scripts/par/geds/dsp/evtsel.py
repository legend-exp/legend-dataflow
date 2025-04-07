import argparse
import json
import time
import warnings
from bisect import bisect_left
from pathlib import Path

import lgdo
import lgdo.lh5 as lh5
import numpy as np
import pygama.math.histogram as pgh
import pygama.pargen.energy_cal as pgc
from dbetto import TextDB
from dbetto.catalog import Props
from pygama.pargen.data_cleaning import generate_cuts, get_keys
from pygama.pargen.dsp_optimize import run_one_dsp

from .....log import build_log
from ....pulser_removal import get_pulser_mask

warnings.filterwarnings(action="ignore", category=RuntimeWarning)


def get_out_data(
    raw_data,
    dsp_data,
    cut_dict,
    e_lower_lim,
    e_upper_lim,
    ecal_pars,
    raw_dict,
    peak,
    final_cut_field="is_valid_cal",
    energy_param="trapTmax",
):
    for outname, info in cut_dict.items():
        outcol = dsp_data.eval(info["expression"], info.get("parameters", None))
        dsp_data.add_column(outname, outcol)

    for outname, info in raw_dict.items():
        outcol = raw_data.eval(info["expression"], info.get("parameters", None))
        raw_data.add_column(outname, outcol)

    final_mask = (
        (dsp_data[energy_param].nda > e_lower_lim)
        & (dsp_data[energy_param].nda < e_upper_lim)
        & (dsp_data[final_cut_field].nda)
    )

    wavefrom_windowed = lgdo.WaveformTable(
        t0=raw_data["waveform_windowed"]["t0"].nda[final_mask],
        t0_units=raw_data["waveform_windowed"]["t0"].attrs["units"],
        dt=raw_data["waveform_windowed"]["dt"].nda[final_mask],
        dt_units=raw_data["waveform_windowed"]["dt"].attrs["units"],
        values=raw_data["waveform_windowed"]["values"].nda[final_mask],
    )
    wavefrom_presummed = lgdo.WaveformTable(
        t0=raw_data["waveform_presummed"]["t0"].nda[final_mask],
        t0_units=raw_data["waveform_presummed"]["t0"].attrs["units"],
        dt=raw_data["waveform_presummed"]["dt"].nda[final_mask],
        dt_units=raw_data["waveform_presummed"]["dt"].attrs["units"],
        values=raw_data["waveform_presummed"]["values"].nda[final_mask],
    )

    out_tbl = lgdo.Table(
        col_dict={
            "waveform_presummed": wavefrom_presummed,
            "waveform_windowed": wavefrom_windowed,
            "presum_rate": lgdo.Array(raw_data["presum_rate"].nda[final_mask]),
            "timestamp": lgdo.Array(raw_data["timestamp"].nda[final_mask]),
            "baseline": lgdo.Array(raw_data["baseline"].nda[final_mask]),
            "daqenergy": lgdo.Array(raw_data["daqenergy"].nda[final_mask]),
            "daqenergy_cal": lgdo.Array(raw_data["daqenergy_cal"].nda[final_mask]),
            "trapTmax_cal": lgdo.Array(
                dsp_data["trapTmax"].nda[final_mask] * ecal_pars
            ),
            "peak": lgdo.Array(np.full(len(np.where(final_mask)[0]), int(peak))),
        }
    )
    return out_tbl, len(np.where(final_mask)[0])


def par_geds_dsp_evtsel() -> None:
    argparser = argparse.ArgumentParser()
    argparser.add_argument("--raw-filelist", help="raw_filelist", type=str)
    argparser.add_argument(
        "--pulser-file", help="pulser-file", type=str, required=False
    )
    argparser.add_argument(
        "-p", "--no-pulse", help="no pulser present", action="store_true"
    )

    argparser.add_argument("--decay-const", help="decay_const", type=str, required=True)
    argparser.add_argument(
        "--raw-cal-curve",
        help="raw calibration curve file(s)",
        type=str,
        nargs="*",
        required=True,
    )

    argparser.add_argument("--log", help="log_file", type=str)
    argparser.add_argument("--configs", help="configs", type=str, required=True)

    argparser.add_argument("--datatype", help="Datatype", type=str, required=True)
    argparser.add_argument("--timestamp", help="Timestamp", type=str, required=True)
    argparser.add_argument("--channel", help="Channel", type=str, required=True)
    argparser.add_argument(
        "--raw-table-name", help="raw table name", type=str, required=True
    )

    argparser.add_argument("--peak-file", help="peak_file", type=str, required=True)
    args = argparser.parse_args()

    configs = TextDB(args.configs, lazy=True).on(args.timestamp, system=args.datatype)
    config_dict = configs["snakemake_rules"]["pars_dsp_peak_selection"]

    log = build_log(config_dict, args.log)

    sto = lh5.LH5Store()
    t0 = time.time()

    dsp_config = config_dict["inputs"]["processing_chain"][args.channel]
    peak_json = config_dict["inputs"]["peak_config"][args.channel]

    peak_dict = Props.read_from(peak_json)
    db_dict = Props.read_from(args.decay_const)

    Path(args.peak_file).parent.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng()
    rand_num = f"{rng.integers(0, 99999):05d}"
    temp_output = f"{args.peak_file}.{rand_num}"
    if peak_dict.pop("run_selection") is True:
        log.debug("Starting peak selection")

        with Path(args.raw_filelist).open() as f:
            files = f.read().splitlines()
        raw_files = sorted(files)

        raw_dict = Props.read_from(args.raw_cal_curve)[args.channel]["pars"][
            "operations"
        ]

        peaks_kev = peak_dict["peaks"]
        kev_widths = peak_dict["kev_widths"]
        cut_parameters = peak_dict["cut_parameters"]
        n_events = peak_dict["n_events"]
        final_cut_field = peak_dict["final_cut_field"]
        energy_parameter = peak_dict.get("energy_parameter", "trapTmax")

        lh5_path = args.raw_table_name

        if not isinstance(kev_widths, list):
            kev_widths = [kev_widths]

        if lh5_path[-1] != "/":
            lh5_path += "/"

        tb = lh5.read(
            lh5_path, raw_files, field_mask=["daqenergy", "t_sat_lo", "timestamp"]
        )

        if args.no_pulse is False:
            mask = get_pulser_mask(
                args.pulser_file,
            )
        else:
            mask = np.full(len(tb), False)

        discharges = tb["t_sat_lo"].nda > 0
        discharge_timestamps = np.where(tb["timestamp"].nda[discharges])[0]
        is_recovering = np.full(len(tb), False, dtype=bool)
        for tstamp in discharge_timestamps:
            is_recovering = is_recovering | np.where(
                (
                    ((tb["timestamp"].nda - tstamp) < 0.01)
                    & ((tb["timestamp"].nda - tstamp) > 0)
                ),
                True,
                False,
            )

        for outname, info in raw_dict.items():
            outcol = tb.eval(info["expression"], info.get("parameters", None))
            tb.add_column(outname, outcol)

        rough_energy = tb["daqenergy_cal"].nda

        masks = {}
        for peak, kev_width in zip(peaks_kev, kev_widths):
            e_mask = (
                (rough_energy > peak - 1.1 * kev_width[0])
                & (rough_energy < peak + 1.1 * kev_width[0])
                & (~mask)
            )
            masks[peak] = np.where(e_mask & (~is_recovering))[0]
            log.debug(f"{len(masks[peak])} events found in energy range for {peak}")

        input_data = lh5.read(
            f"{lh5_path}", raw_files, n_rows=10000, idx=np.where(~mask)[0]
        )

        if isinstance(dsp_config, str):
            dsp_config = Props.read_from(dsp_config)

        dsp_config["outputs"] = [
            *get_keys(dsp_config["outputs"], cut_parameters),
            energy_parameter,
        ]

        log.debug("Processing data")
        tb_data = run_one_dsp(input_data, dsp_config, db_dict=db_dict)

        if cut_parameters is not None:
            cut_dict = generate_cuts(tb_data, cut_parameters)
            log.debug(f"Cuts are calculated: {json.dumps(cut_dict, indent=2)}")
        else:
            cut_dict = None

        pk_dicts = {}
        for peak, kev_width in zip(peaks_kev, kev_widths):
            pk_dicts[peak] = {
                "idxs": (masks[peak],),
                "n_rows_read": 0,
                "obj_buf_start": 0,
                "obj_buf": None,
                "kev_width": kev_width,
            }

        for file in raw_files:
            log.debug(Path(file).name)
            for peak, peak_dict in pk_dicts.items():
                if peak_dict["idxs"] is not None:
                    # idx is a long continuous array
                    n_rows_i = sto.read_n_rows(lh5_path, file)
                    # find the length of the subset of idx that contains indices
                    # that are less than n_rows_i
                    n_rows_to_read_i = bisect_left(peak_dict["idxs"][0], n_rows_i)
                    # now split idx into idx_i and the remainder
                    idx_i = (peak_dict["idxs"][0][:n_rows_to_read_i],)
                    peak_dict["idxs"] = (
                        peak_dict["idxs"][0][n_rows_to_read_i:] - n_rows_i,
                    )
                    if len(idx_i[0]) > 0:
                        peak_dict["obj_buf"] = lh5.read(
                            lh5_path,
                            file,
                            start_row=0,
                            idx=idx_i,
                            obj_buf=peak_dict["obj_buf"],
                            obj_buf_start=peak_dict["obj_buf_start"],
                        )
                        n_rows_read_i = len(peak_dict["obj_buf"])

                        peak_dict["n_rows_read"] += n_rows_read_i
                        log.debug(f"{peak}: {peak_dict['n_rows_read']}")
                        peak_dict["obj_buf_start"] += n_rows_read_i
                    if peak_dict["n_rows_read"] >= 10000 or file == raw_files[-1]:
                        if "e_lower_lim" not in peak_dict:
                            tb_out = run_one_dsp(
                                peak_dict["obj_buf"], dsp_config, db_dict=db_dict
                            )
                            energy = tb_out[energy_parameter].nda

                            init_bin_width = (
                                2
                                * (
                                    np.nanpercentile(energy, 75)
                                    - np.nanpercentile(energy, 25)
                                )
                                * len(energy) ** (-1 / 3)
                            )

                            if init_bin_width > 2:
                                init_bin_width = 2

                            hist, bins, var = pgh.get_hist(
                                energy,
                                range=(
                                    np.floor(np.nanpercentile(energy, 1)),
                                    np.ceil(np.nanpercentile(energy, 99)),
                                ),
                                dx=init_bin_width,
                            )
                            peak_loc = pgh.get_bin_centers(bins)[np.nanargmax(hist)]

                            peak_top_pars = pgc.hpge_fit_energy_peak_tops(
                                hist,
                                bins,
                                var,
                                [peak_loc],
                                n_to_fit=7,
                            )[0][0]
                            try:
                                mu = peak_top_pars[0]
                                if mu > np.nanmax(bins) or mu < np.nanmin(bins):
                                    raise ValueError
                            except Exception:
                                mu = np.nan
                            if mu is None or np.isnan(mu):
                                log.debug("Fit failed, using max guess")
                                rough_adc_to_kev = peak / peak_loc
                                e_lower_lim = (
                                    peak_loc
                                    - (1.5 * peak_dict["kev_width"][0])
                                    / rough_adc_to_kev
                                )
                                e_upper_lim = (
                                    peak_loc
                                    + (1.5 * peak_dict["kev_width"][1])
                                    / rough_adc_to_kev
                                )
                                hist, bins, var = pgh.get_hist(
                                    energy,
                                    range=(int(e_lower_lim), int(e_upper_lim)),
                                    dx=init_bin_width,
                                )
                                mu = pgh.get_bin_centers(bins)[np.nanargmax(hist)]

                            updated_adc_to_kev = peak / mu
                            e_lower_lim = (
                                mu - (peak_dict["kev_width"][0]) / updated_adc_to_kev
                            )
                            e_upper_lim = (
                                mu + (peak_dict["kev_width"][1]) / updated_adc_to_kev
                            )
                            log.info(
                                f"{peak}: lower lim is :{e_lower_lim}, upper lim is {e_upper_lim}"
                            )
                            peak_dict["e_lower_lim"] = e_lower_lim
                            peak_dict["e_upper_lim"] = e_upper_lim
                            peak_dict["ecal_par"] = updated_adc_to_kev

                            out_tbl, n_wfs = get_out_data(
                                peak_dict["obj_buf"],
                                tb_out,
                                cut_dict,
                                e_lower_lim,
                                e_upper_lim,
                                peak_dict["ecal_par"],
                                raw_dict,
                                int(peak),
                                final_cut_field=final_cut_field,
                                energy_param=energy_parameter,
                            )
                            lh5.write(
                                out_tbl,
                                name=lh5_path,
                                lh5_file=temp_output,
                                wo_mode="a",
                            )
                            peak_dict["obj_buf"] = None
                            peak_dict["obj_buf_start"] = 0
                            peak_dict["n_events"] = n_wfs
                            log.debug(
                                f"found {peak_dict['n_events']} events for {peak}"
                            )
                        else:
                            if (
                                peak_dict["obj_buf"] is not None
                                and len(peak_dict["obj_buf"]) > 0
                            ):
                                tb_out = run_one_dsp(
                                    peak_dict["obj_buf"], dsp_config, db_dict=db_dict
                                )
                                out_tbl, n_wfs = get_out_data(
                                    peak_dict["obj_buf"],
                                    tb_out,
                                    cut_dict,
                                    peak_dict["e_lower_lim"],
                                    peak_dict["e_upper_lim"],
                                    peak_dict["ecal_par"],
                                    raw_dict,
                                    int(peak),
                                    final_cut_field=final_cut_field,
                                    energy_param=energy_parameter,
                                )
                                peak_dict["n_events"] += n_wfs
                                lh5.write(
                                    out_tbl,
                                    name=lh5_path,
                                    lh5_file=temp_output,
                                    wo_mode="a",
                                )
                                peak_dict["obj_buf"] = None
                                peak_dict["obj_buf_start"] = 0
                                log.debug(
                                    f"found {peak_dict['n_events']} events for {peak}"
                                )
                                if peak_dict["n_events"] >= n_events:
                                    peak_dict["idxs"] = None
                                    log.debug(
                                        f"{peak} has reached the required number of events"
                                    )

    else:
        Path(temp_output).touch()

    log.debug(f"event selection completed in {time.time() - t0} seconds")
    Path(temp_output).rename(args.peak_file)
