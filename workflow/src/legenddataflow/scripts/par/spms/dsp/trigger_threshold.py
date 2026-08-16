"""Console scripts ``par-spms-dsp-trg-thr`` and ``-multi``: determine SiPM
trigger thresholds from the noise distribution of DSP-processed raw
waveforms."""

from __future__ import annotations

import argparse

import hist
import lh5
import numpy as np
from dbetto import AttrsDict, Props, utils
from dspeed import build_dsp
from legenddataflowscripts.utils import (
    build_log,
    cfgtools,
    check_input_files,
    get_rule_config,
    prepare_output_paths,
)


def get_channel_trg_thr(df_configs, sipm_name, dsp_db, raw_file, raw_table_name, log):
    log.debug("reading in the configuration files")
    config = df_configs.inputs
    dsp_config = utils.load_dict(
        cfgtools.get_channel_config(config.processing_chain, sipm_name)
    )
    settings = AttrsDict(
        utils.load_dict(cfgtools.get_channel_config(config.settings, sipm_name))
    )

    # get DSP database from overrides
    _db_dict = Props.read_from(dsp_db).get(sipm_name, {})

    fwhm = None

    # read raw file list
    log.debug("reading in the raw waveforms")
    if len(lh5.ls(raw_file, f"{raw_table_name}/waveform_bit_drop")) == 0:
        msg = (
            f"could not find waveform '{raw_table_name}/waveform_bit_drop' "
            f"in {raw_file}, returning null pars"
        )
        log.warning(msg)

    else:
        data = lh5.read(
            raw_table_name,
            raw_file,
            field_mask=["waveform_bit_drop"],
            n_rows=settings.n_events,
        )

        if len(data) == 0:
            msg = (
                f"could not find any waveforms '{raw_table_name}/waveform_bit_drop' "
                f"in {raw_file}, returning null pars"
            )
            log.warning(msg)

        elif len(data) < settings.n_events:
            msg = (
                f"found only {len(data)} waveforms in "
                f"'{raw_table_name}/waveform_bit_drop' of {raw_file}, need at "
                f"least {settings.n_events} to build the histogram"
            )
            raise RuntimeError(msg)
        else:
            # run the DSP with the provided configuration
            log.debug("running the DSP chain")
            dsp_output = build_dsp(
                raw_in=data, dsp_config=dsp_config, database=_db_dict
            )

            log.debug("analyzing DSP outputs")
            # get output of the current processor
            wf_current = dsp_output.wf_current.values.view_as("np").flatten()
            # determine a cutoff for the histogram used to extract the FWHM
            low_cutoff, high_cutoff = np.quantile(wf_current, [0.005, 0.995])

            # determine hist edges with Friedmann Diaconis Estimator
            bin_edges = np.histogram_bin_edges(
                wf_current, bins="fd", range=(low_cutoff, high_cutoff)
            )

            # make histogram of the curr values
            h = hist.new.Variable(bin_edges).Double().fill(wf_current)

            # determine FWHM
            counts = h.view()
            idx_over_half = np.where(counts >= np.max(counts) / 2)[0]

            edges = h.axes[0].edges
            fwhm = edges[idx_over_half[-1]] - edges[idx_over_half[0]]

            if fwhm <= 0:
                msg = (
                    "determined FWHM of the baseline derivative distribution "
                    f"is not positive: {fwhm:.3f}"
                )
                raise RuntimeError(msg)

            return fwhm
    return None


def par_spms_dsp_trg_thr() -> None:
    # CLI interface
    argparser = argparse.ArgumentParser()
    argparser.add_argument("--raw-file", required=True)
    argparser.add_argument("--raw-table-name", required=True)
    argparser.add_argument("--output-file", required=True)
    argparser.add_argument("--config-path", required=True)
    argparser.add_argument("--datatype", required=True)
    argparser.add_argument("--timestamp", required=True)
    argparser.add_argument("--sipm-name", required=True)
    argparser.add_argument("--dsp-db", nargs="*", default=[])
    argparser.add_argument("--logfile")
    args = argparser.parse_args()

    # dataflow configs
    df_configs = get_rule_config(
        args.config_path, "pars_spms_dsp_trg_thr", args.timestamp, args.datatype
    )

    # setup logging
    log = build_log(df_configs, args.logfile)

    check_input_files(args.raw_file, "--raw-file")
    check_input_files(args.dsp_db, "--dsp-db")
    prepare_output_paths(args.output_file)

    fwhm = get_channel_trg_thr(
        df_configs,
        args.sipm_name,
        args.dsp_db,
        args.raw_file,
        args.raw_table_name,
        log,
    )
    msg = f"writing out baseline_curr_fwhm = {fwhm}"
    log.debug(msg)
    Props.write_to(
        args.output_file,
        {"baseline_curr_fwhm": float(fwhm) if fwhm is not None else fwhm},
    )


def par_spms_dsp_trg_thr_multi() -> None:
    # CLI interface
    argparser = argparse.ArgumentParser()
    argparser.add_argument("--raw-file", required=True)
    argparser.add_argument("--raw-table-names", required=True, nargs="*")
    argparser.add_argument("--output-file", required=True)
    argparser.add_argument("--config-path", required=True)
    argparser.add_argument("--datatype", required=True)
    argparser.add_argument("--timestamp", required=True)
    argparser.add_argument("--sipm-names", required=True, nargs="*")
    argparser.add_argument("--dsp-db", nargs="*", default=[])
    argparser.add_argument("--logfile")
    args = argparser.parse_args()

    # dataflow configs
    df_configs = get_rule_config(
        args.config_path, "pars_spms_dsp_trg_thr", args.timestamp, args.datatype
    )

    # setup logging
    log = build_log(df_configs, args.logfile)

    check_input_files(args.raw_file, "--raw-file")
    check_input_files(args.dsp_db, "--dsp-db")
    prepare_output_paths(args.output_file)

    if len(args.sipm_names) != len(args.raw_table_names):
        msg = (
            f"--sipm-names ({len(args.sipm_names)} entries) and "
            f"--raw-table-names ({len(args.raw_table_names)} entries) must "
            "have the same length"
        )
        raise ValueError(msg)

    out_dict = {}
    for sipm_name, raw_table_name in zip(
        args.sipm_names, args.raw_table_names, strict=True
    ):
        fwhm = get_channel_trg_thr(
            df_configs,
            sipm_name,
            args.dsp_db,
            args.raw_file,
            raw_table_name,
            log,
        )

        msg = f"baseline_curr_fwhm for {sipm_name} = {fwhm}"
        log.debug(msg)
        out_dict[sipm_name] = {
            "baseline_curr_fwhm": float(fwhm) if fwhm is not None else fwhm
        }
    Props.write_to(
        args.output_file,
        out_dict,
    )
