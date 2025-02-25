import argparse
from pathlib import Path

import hist
import numpy as np
from dbetto import AttrsDict, Props, TextDB, utils
from dspeed import build_processing_chain
from lgdo import lh5

from ..... import cfgtools
from .....log import build_log


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
    df_configs = (
        TextDB(args.config_path, lazy=True)
        .on(args.timestamp, system=args.datatype)
        .snakemake_rules.pars_spms_dsp_trg_thr
    )

    # setup logging
    log = build_log(df_configs, args.logfile)

    log.debug("reading in the configuration files")
    config = df_configs.inputs
    dsp_config = utils.load_dict(
        cfgtools.get_channel_config(config.processing_chain, args.sipm_name)
    )
    settings = AttrsDict(
        utils.load_dict(cfgtools.get_channel_config(config.settings, args.sipm_name))
    )

    # read raw file list
    log.debug("reading in the raw waveforms")
    data = lh5.read(
        args.raw_table_name,
        args.raw_file,
        field_mask=["waveform_bit_drop"],
        n_rows=settings.n_events,
    )

    # get DSP database from overrides
    _db_dict = Props.read_from(args.dsp_db).get(args.sipm_name, {})

    # run the DSP with the provided configuration
    log.debug("running the DSP chain")
    chain, _, dsp_output = build_processing_chain(data, dsp_config, db_dict=_db_dict)
    chain.execute()

    log.debug("analyzing DSP outputs")
    # get output of the current processor
    wf_current = dsp_output.wf_current.values.view_as("np").flatten()
    # determine a cutoff for the histogram used to extract the FWHM
    low_cutoff, high_cutoff = np.quantile(wf_current, [0.005, 0.995])

    # make histogram of the curr values
    h = (
        hist.new.Regular(settings.n_baseline_bins, low_cutoff, high_cutoff)
        .Double()
        .fill(wf_current)
    )

    # determine FWHM
    counts = h.view()
    idx_over_half = np.where(counts >= np.max(counts) / 2)[0]

    edges = h.axes[0].edges
    fwhm = edges[idx_over_half[-1]] - edges[idx_over_half[0]]

    if fwhm <= 0:
        msg = "determined FWHM of baseline derivative distribution is zero or negative"
        raise RuntimeError(msg)

    log.debug(f"writing out baseline_curr_fwhm = {fwhm}")
    Path(args.output_file).parent.mkdir(parents=True, exist_ok=True)
    Props.write_to(args.output_file, {"baseline_curr_fwhm": float(fwhm)})
