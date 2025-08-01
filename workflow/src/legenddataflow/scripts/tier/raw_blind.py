"""
This script takes in raw data, applies the calibration to the daqenergy
and uses this to blind the data in a window of Qbb +- 25 keV. It copies over all
channels in a raw file, removing those events that fall within the ROI for Ge detectors
that have a daqenergy calibration curve and are not anti-coincidence only (AC). It removes
the whole event from all of the Ge and SiPM channels.

In the Snakemake dataflow, this script only runs if the checkfile is found on disk,
but this is controlled by the Snakemake flow (presumably an error is thrown if the file
is not found). This script itself does not check for the existence of such a file.
"""

from __future__ import annotations

import argparse
import time

import hdf5plugin
import numexpr as ne
import numpy as np
from dbetto.catalog import Props
from legenddataflowscripts.utils import alias_table, build_log
from legendmeta import LegendMetadata, TextDB
from lgdo import lh5

filter_map = {
    "zstd": hdf5plugin.Zstd(),
    "blosc": hdf5plugin.Blosc(),
    "lz4": hdf5plugin.LZ4(),
    "bitshuffle": hdf5plugin.Bitshuffle(),
    "bzip2": hdf5plugin.BZip2(),
    # Add other filters from hdf5plugin if needed
}


def build_tier_raw_blind() -> None:
    argparser = argparse.ArgumentParser()
    argparser.add_argument("--input", help="input file", type=str)
    argparser.add_argument("--output", help="output file", type=str)
    argparser.add_argument(
        "--blind-curve", help="blinding curves file", type=str, required=True, nargs="*"
    )
    argparser.add_argument("--datatype", help="Datatype", type=str, required=True)
    argparser.add_argument("--timestamp", help="Timestamp", type=str, required=True)
    argparser.add_argument("--configs", help="config file", type=str)
    argparser.add_argument("--chan-maps", help="chan map", type=str)
    argparser.add_argument("--metadata", help="metadata", type=str)
    argparser.add_argument("--log", help="log file", type=str)
    argparser.add_argument("--alias-table", help="Alias table", type=str, default=None)
    args = argparser.parse_args()

    configs = TextDB(args.configs, lazy=True)
    config_dict = configs.on(args.timestamp, system=args.datatype)["snakemake_rules"][
        "tier_raw_blind"
    ]

    log = build_log(config_dict, args.log)

    hdf_settings = Props.read_from(config_dict["settings"])["hdf5_settings"]

    if "compression" in hdf_settings:
        compression = hdf_settings["compression"]
        if compression in filter_map:
            hdf_settings["compression"] = filter_map[compression]

    blinding_settings = Props.read_from(config_dict["config"])

    centroid = blinding_settings["centroid_in_keV"]  # keV
    width = blinding_settings["width_in_keV"]  # keV

    # list of all channels and objects in the raw file
    all_channels = lh5.ls(args.input)

    # list of Ge channels and SiPM channels with associated metadata
    legendmetadata = LegendMetadata(args.metadata, lazy=True)
    chmap = legendmetadata.channelmap(args.timestamp)
    chans = {
        system: chan_dict.map("daq.rawid")
        for system, chan_dict in chmap.map("system", unique=False).items()
    }

    main_channels = (
        list(chans["geds"])
        + list(chans["spms"])
        + list(chans["auxs"])
        + list(chans["blsns"])
        + list(chans["puls"])
    )

    # rows that need blinding
    toblind = np.array([])

    log.info("Blinding Ge channels")
    start = time.time()

    # first, loop through the Ge detector channels, calibrate them and look for events that should be blinded
    for chnum in chans["geds"]:
        # skip Ge detectors that are anti-coincidence only or not able to be blinded for some other reason
        if chans["geds"][chnum]["analysis"]["is_blinded"] is False:
            continue

        # load in just the daqenergy for now
        daqenergy = lh5.read(f"ch{chnum}/raw/daqenergy", args.input)

        # read in calibration curve for this channel
        blind_curve = Props.read_from(args.blind_curve)[chmap.map("daq.rawid").name][
            "pars"
        ]["operations"]

        # calibrate daq energy using pre existing curve
        daqenergy_cal = ne.evaluate(
            blind_curve["daqenergy_cal"]["expression"],
            local_dict=dict(
                daqenergy=daqenergy, **blind_curve["daqenergy_cal"]["parameters"]
            ),
        )

        # figure out which event indices should be blinded
        toblind = np.append(
            toblind,
            np.nonzero(np.abs(np.asarray(daqenergy_cal) - centroid) <= width)[0],
        )

    # remove duplicates
    toblind = np.unique(toblind)

    # total number of events (from last Ge channel loaded, should be same for all Ge channels)
    allind = np.arange(len(daqenergy))

    # gets events that should not be blinded
    tokeep = allind[np.logical_not(np.isin(allind, toblind))]

    for channel in all_channels:
        try:
            chnum = int(channel[2::])
        except ValueError:
            # if this isn't an interesting channel, just copy it to the output file
            chobj = lh5.read(channel, args.input, decompress=False)
            lh5.write(
                chobj,
                channel,
                lh5_file=args.output,
                wo_mode="w",
                **hdf_settings,
            )
            continue

        if chnum not in main_channels:
            # if this is a PMT or not included for some reason, just copy it to the output file
            chobj = lh5.read(channel + "/raw", args.input, decompress=False)
            lh5.write(
                chobj,
                group=channel,
                name="raw",
                lh5_file=args.output,
                wo_mode="w",
                **hdf_settings,
            )
            continue

        # the rest should be the Ge and SiPM channels that need to be blinded

        # read in all of the data but only for the unblinded events
        blinded_chobj = lh5.read(
            channel + "/raw", args.input, idx=tokeep, decompress=False
        )

        # now write the blinded data for this channel
        lh5.write(
            blinded_chobj,
            group=channel,
            name="raw",
            lh5_file=args.output,
            wo_mode="w",
            **hdf_settings,
        )

    log.info("Finished blinding Ge channels")
    msg = f"Time taken: {time.time() - start:.2f} seconds"
    log.info(msg)
    if args.alias_table is not None:
        log.info("Creating alias table")
        alias_table(args.output, args.alias_table)
