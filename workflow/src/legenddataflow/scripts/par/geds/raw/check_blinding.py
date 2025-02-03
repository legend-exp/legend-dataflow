"""
This script checks that the blinding for a particular channel is still valid,
it does this by taking the calibration curve stored in the overrides, applying it
to the daqenergy, running a peak search over the calibrated energy and checking that
there are peaks within 5keV of the 583 and 2614 peaks. If the detector is in ac mode
then it will skip the check.
"""

import argparse
import pickle as pkl
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numexpr as ne
import numpy as np
from dbetto import TextDB
from dbetto.catalog import Props
from legendmeta import LegendMetadata
from lgdo import lh5
from pygama.math.histogram import get_hist
from pygama.pargen.energy_cal import get_i_local_maxima

from ....log import build_log

mpl.use("Agg")


def par_geds_raw_blindcheck() -> None:
    argparser = argparse.ArgumentParser()
    argparser.add_argument("--files", help="files", nargs="*", type=str)
    argparser.add_argument("--output", help="output file", type=str)
    argparser.add_argument("--plot_file", help="plot file", type=str)
    argparser.add_argument(
        "--blind_curve", help="blinding curves file", nargs="*", type=str
    )
    argparser.add_argument("--datatype", help="Datatype", type=str, required=True)
    argparser.add_argument("--timestamp", help="Timestamp", type=str, required=True)
    argparser.add_argument("--configs", help="config file", type=str)
    argparser.add_argument("--channel", help="channel", type=str)
    argparser.add_argument("--metadata", help="channel", type=str)
    argparser.add_argument("--log", help="log file", type=str)
    args = argparser.parse_args()

    configs = TextDB(args.configs, lazy=True).on(args.timestamp, system=args.datatype)
    config_dict = configs["snakemake_rules"]["tier_raw_blindcheck"]

    log = build_log(config_dict, args.log)

    # get the usability status for this channel
    chmap = (
        LegendMetadata(args.metadata, lazy=True)
        .channelmap(args.timestamp)
        .map("daq.rawid")
    )
    det_status = chmap[int(args.channel[2:])]["analysis"]["is_blinded"]

    # read in calibration curve for this channel
    blind_curve = Props.read_from(args.blind_curve)[args.channel]["pars"]["operations"]

    # load in the data
    daqenergy = lh5.read(f"{args.channel}/raw/daqenergy", sorted(args.files))[
        0
    ].view_as("np")

    # calibrate daq energy using pre existing curve
    daqenergy_cal = ne.evaluate(
        blind_curve["daqenergy_cal"]["expression"],
        local_dict=dict(
            daqenergy=daqenergy, **blind_curve["daqenergy_cal"]["parameters"]
        ),
    )

    # bin with 1 keV bins and get maxs
    hist, bins, var = get_hist(daqenergy_cal, np.arange(0, 3000, 1))
    maxs = get_i_local_maxima(hist, delta=25)
    log.info(f"peaks found at : {maxs}")

    # plot the energy spectrum to check calibration
    fig = plt.figure(figsize=(8, 10))
    ax = plt.subplot(211)
    ax.hist(daqenergy_cal, bins=np.arange(0, 3000, 1), histtype="step")
    ax.set_ylabel("counts")
    ax.set_yscale("log")
    ax2 = plt.subplot(212)
    ax2.hist(
        daqenergy_cal,
        bins=np.arange(2600, 2630, 1 * blind_curve["daqenergy_cal"]["parameters"]["a"]),
        histtype="step",
    )
    ax2.set_xlabel("energy (keV)")
    ax2.set_ylabel("counts")
    plt.suptitle(args.channel)
    with Path(args.plot_file).open("wb") as w:
        pkl.dump(fig, w, protocol=pkl.HIGHEST_PROTOCOL)
    plt.close()

    # check for peaks within +- 5keV of  2614 and 583 to ensure blinding still
    # valid and if so create file else raise error.  if detector is in ac mode it
    # will always pass this check
    if (
        np.any(np.abs(maxs - 2614) < 5) and np.any(np.abs(maxs - 583) < 5)
    ) or det_status is False:
        Path(args.output).mkdir(parents=True, exist_ok=True)
        Props.write_to(
            args.output,
            {
                "threshold_adc": np.nanmin(daqenergy),
                "threshold_kev": np.nanmin(daqenergy_cal),
            },
        )
    else:
        msg = "peaks not found in daqenergy"
        raise RuntimeError(msg)
