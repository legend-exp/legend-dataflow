from __future__ import annotations

import argparse
import time
from pathlib import Path

import numpy as np
import pygama.math.distributions as pmd  # noqa: F401
from dbetto.catalog import Props
from legenddataflowscripts.utils import build_log
from lgdo import Array, Table, WaveformTable, lh5
from pygama.evt.build_tcm import _concat_tables
from pygama.pargen.data_cleaning import generate_cuts
from pygama.pargen.dplms_ge_dict import dplms_ge_dict
from pygama.pargen.dsp_optimize import run_one_dsp

from legenddataflow.methods.FileKey import ChannelProcKey
from legenddataflow.scripts.par.geds.pht.util import (
    get_run_dict,
    save_dict_to_files,
    split_files_by_run,
)


def filter_table(raw_data_tables, final_masks):

    waveform_windowed = WaveformTable(
        t0=np.concatenate(
            [
                tbl["waveform_windowed"]["t0"].nda[mask]
                for tbl, mask in zip(raw_data_tables, final_masks, strict=True)
            ]
        ),
        t0_units=raw_data_tables[0]["waveform_windowed"]["t0"].attrs["units"],
        dt=np.concatenate(
            [
                tbl["waveform_windowed"]["dt"].nda[mask]
                for tbl, mask in zip(raw_data_tables, final_masks, strict=True)
            ]
        ),
        dt_units=raw_data_tables[0]["waveform_windowed"]["dt"].attrs["units"],
        values=np.concatenate(
            [
                tbl["waveform_windowed"]["values"].nda[mask]
                for tbl, mask in zip(raw_data_tables, final_masks, strict=True)
            ]
        ),
    )
    waveform_presummed = WaveformTable(
        t0=np.concatenate(
            [
                tbl["waveform_presummed"]["t0"].nda[mask]
                for tbl, mask in zip(raw_data_tables, final_masks, strict=True)
            ]
        ),
        t0_units=raw_data_tables[0]["waveform_presummed"]["t0"].attrs["units"],
        dt=np.concatenate(
            [
                tbl["waveform_presummed"]["dt"].nda[mask]
                for tbl, mask in zip(raw_data_tables, final_masks, strict=True)
            ]
        ),
        dt_units=raw_data_tables[0]["waveform_presummed"]["dt"].attrs["units"],
        values=np.concatenate(
            [
                tbl["waveform_presummed"]["values"].nda[mask]
                for tbl, mask in zip(raw_data_tables, final_masks, strict=True)
            ]
        ),
    )

    col_dict = {
        "waveform_presummed": waveform_presummed,
        "waveform_windowed": waveform_windowed,
        "presum_rate": Array(
            np.concatenate(
                [
                    tbl["presum_rate"].nda[mask]
                    for tbl, mask in zip(raw_data_tables, final_masks, strict=True)
                ]
            )
        ),
        "timestamp": Array(
            np.concatenate(
                [
                    tbl["timestamp"].nda[mask]
                    for tbl, mask in zip(raw_data_tables, final_masks, strict=True)
                ]
            )
        ),
        "baseline": Array(
            np.concatenate(
                [
                    tbl["baseline"].nda[mask]
                    for tbl, mask in zip(raw_data_tables, final_masks, strict=True)
                ]
            )
        ),
        "daqenergy": Array(
            np.concatenate(
                [
                    tbl["daqenergy"].nda[mask]
                    for tbl, mask in zip(raw_data_tables, final_masks, strict=True)
                ]
            )
        ),
    }
    return Table(col_dict)


def par_geds_psp_dplms() -> None:
    """Compute DPLMS optimal filter coefficients for HPGe detectors.
    This is the partition version of dplms, it combines multiple runs
    with an equal number of fft and cal events from each.

    CLI entry point registered as ``par-geds-psp-dplms``.  The *Discrete
    Prolate Laplacian Maximum Signal* (DPLMS) filter is constructed from two
    data samples:

    1. **FFT baselines** - low-energy (``daqenergy <= 10`` ADC) events from
       dedicated FFT run files used to characterise the noise power spectrum.
    2. **Calibration peak events** - gamma-line events read from the
       *peak-file* produced by :func:`par_geds_dsp_evtsel`.

    The coefficients are computed by
    :func:`pygama.pargen.dplms_ge_dict.dplms_ge_dict` and written to an LH5
    file at *lh5-path* (so that the DSP chain can load them with
    ``loadlh5(...)``).  The DSP parameter database is updated with the filter
    configuration and written to *dsp-pars*.

    Notes
    -----
    **Command-line arguments**

    ``--fft-raw-filelist`` : str
        Path to a text file listing the FFT raw LH5 input files.
    ``--peak-file`` : str
        LH5 file containing preselected calibration peak events.
    ``--inplots`` : str, optional
        Existing pickle plot file to update with DPLMS plots.
    ``--database`` : str
        Path to the existing DSP parameter database (JSON/YAML).
    ``--log`` : str, optional
        Path to the log file.
    ``--log-config`` : str, optional
        Logging configuration file.
    ``--processing-chain`` : list of str
        Processing chain configuration file(s).
    ``--config-file`` : list of str
        DPLMS configuration file(s).  Must contain ``run_dplms`` (bool),
        ``n_baselines`` (int), and ``peaks_kev`` (list).
    ``--channel`` : str
        Channel identifier; used as the HDF5 group name in *lh5-path*.
    ``--raw-table-name`` : str
        LH5 table path within the raw and peak files.
    ``--dsp-pars`` : str
        Output path for the updated DSP parameter database (JSON/YAML).
    ``--lh5-path`` : str
        Output LH5 file path for the DPLMS filter coefficients.
    ``--plot-path`` : str, optional
        Output path for diagnostic plots (pickle).
    """
    argparser = argparse.ArgumentParser()
    argparser.add_argument(
        "--fft-raw-filelists",
        help="fft_raw_filelists",
        type=str,
        nargs="*",
        required=True,
    )
    argparser.add_argument(
        "--peak-files", help="peak files", nargs="*", type=str, required=True
    )
    argparser.add_argument("--inplots", help="in_plot_path", nargs="*", type=str)
    argparser.add_argument(
        "--database", help="database", type=str, nargs="*", required=True
    )

    argparser.add_argument("--log", help="log_file", type=str)
    argparser.add_argument(
        "--log-config", help="Log config file", type=str, required=False, default={}
    )

    argparser.add_argument(
        "--processing-chain",
        help="Processing chain config",
        type=str,
        nargs="*",
        required=True,
    )
    argparser.add_argument(
        "--config-file", help="Config file", type=str, nargs="*", required=True
    )

    argparser.add_argument("--channel", help="channel", type=str, required=True)
    argparser.add_argument(
        "--raw-table-name", help="raw table name", type=str, required=True
    )

    argparser.add_argument(
        "--dsp-pars", help="dsp_pars", type=str, nargs="*", required=True
    )
    argparser.add_argument(
        "--lh5-path", help="lh5_path", type=str, nargs="*", required=True
    )
    argparser.add_argument("--plot-path", help="plot_path", nargs="*", type=str)

    args = argparser.parse_args()

    dsp_config = Props.read_from(args.processing_chain)
    log = build_log(args.log_config, args.log)

    t0 = time.time()

    dplms_dict = Props.read_from(args.config_file)

    db_dicts = get_run_dict(args.database)
    db_dict = db_dicts[next(iter(db_dicts))]

    if dplms_dict["run_dplms"] is True:
        fft_runs, _ = split_files_by_run(args.fft_raw_filelists)

        n_runs = len(fft_runs)
        n_per_run = dplms_dict["n_baselines"] // n_runs

        t0 = time.time()
        log.info("\nLoad fft data")
        fft_tables = []
        eidxs_list = []
        raw_tables = []
        for _, files in fft_runs.items():
            energies = lh5.read_as(
                f"{args.raw_table_name}/daqenergy", files, library="np"
            )
            eidxs = np.where(energies <= 10)[0]
            eidxs_list.append(eidxs)
            raw_fft = lh5.read(
                args.raw_table_name,
                files,
                n_rows=int(dplms_dict["n_baselines"] * 1.5),
                idx=eidxs,
            )
            raw_tables.append(raw_fft)
            dsp_fft = run_one_dsp(raw_fft, dsp_config, db_dict=db_dict)
            fft_tables.append(dsp_fft)

        all_fft = _concat_tables(fft_tables)

        cut_dict = generate_cuts(all_fft, cut_dict=dplms_dict["bls_cut_pars"])
        log.debug("Cuts are %s", cut_dict)

        cut_list = []
        for tbl in fft_tables:
            idxs = np.full(len(tbl), True, dtype=bool)
            for outname, info in cut_dict.items():
                outcol = tbl.eval(info["expression"], info.get("parameters", None))
                tbl.add_column(outname, outcol)
            for cut in cut_dict:
                idxs = tbl[cut].nda & idxs
            cut_list.append(np.where(idxs)[0][:n_per_run])

        raw_fft = filter_table(raw_tables, cut_list)

        log.debug("Applied Cuts")

        t1 = time.time()
        msg = f"Time to load fft data {(t1 - t0):.2f} s, total events {len(raw_fft)}"
        log.info(msg)

        log.info("\nRunning event selection")
        peaks_kev = np.array(dplms_dict["peaks_kev"])
        # kev_widths = [tuple(kev_width) for kev_width in dplms_dict["kev_widths"]]

        peaks_rounded = [int(peak) for peak in peaks_kev]

        n_per_run_signals = dplms_dict["n_signals"] // n_runs
        cal_tbs = []
        for file in args.peak_files:
            peaks = lh5.read_as(f"{args.raw_table_name}/peak", file, library="np")
            ids = np.isin(peaks, peaks_rounded)
            peaks = peaks[ids]

            cal_tbs.append(
                lh5.read(args.raw_table_name, file, idx=ids, n_rows=n_per_run_signals)
            )

        raw_cal = filter_table(cal_tbs, [np.ones(len(t), dtype=bool) for t in cal_tbs])

        msg = f"Time to run event selection {(time.time() - t1):.2f} s, total events {len(raw_cal)}"
        log.info(msg)

        if isinstance(dsp_config, str | list):
            dsp_config = Props.read_from(dsp_config)

        out_dict, plot_dict = dplms_ge_dict(
            raw_fft,
            raw_cal,
            dsp_config,
            db_dict,
            dplms_dict,
            fom_func=eval(dplms_dict.get("fom_func", "pmd.gauss_on_step")),
            display=1 if args.plot_path else 0,
        )

        coeffs = out_dict["dplms"].pop("coefficients")
        dplms_pars = Table(col_dict={"coefficients": Array(coeffs)})
        msg = f"DPLMS creation finished in {(time.time() - t0) / 60} minutes"
        log.info(msg)
    else:
        out_dict = {}
        dplms_pars = Table(col_dict={"coefficients": Array([])})
        plot_dict = {}

    inplot_dict = get_run_dict(args.inplots) if args.inplots else {}

    for _, entry in inplot_dict.items():
        entry.update({"dplms": plot_dict})

    lh5_path_by_ts = {
        ChannelProcKey.get_filekey_from_pattern(Path(p).name).timestamp: p
        for p in args.lh5_path
    }

    for ts, db_dict in db_dicts.items():
        db_dict.update(out_dict)
        if dplms_dict["run_dplms"] is True:
            db_dict.setdefault("dplms", {})["coefficients"] = (
                f"loadlh5('{lh5_path_by_ts[ts]}', '{args.channel}/dplms/coefficients')"
            )

    for outfile in args.lh5_path:
        Path(outfile).parent.mkdir(parents=True, exist_ok=True)
        lh5.write(
            Table(col_dict={"dplms": dplms_pars}),
            name=args.channel,
            lh5_file=outfile,
            wo_mode="overwrite",
        )

    save_dict_to_files(args.dsp_pars, db_dicts)

    if args.plot_path:
        save_dict_to_files(args.plot_path, inplot_dict)
