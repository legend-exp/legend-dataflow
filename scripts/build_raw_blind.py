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

import argparse
import logging
from pathlib import Path

import numexpr as ne
import numpy as np
from legendmeta import LegendMetadata, TextDB
from legendmeta.catalog import Props
from lgdo import lh5

argparser = argparse.ArgumentParser()
argparser.add_argument("--input", help="input file", type=str)
argparser.add_argument("--output", help="output file", type=str)
argparser.add_argument(
    "--blind_curve", help="blinding curves file", type=str, required=True, nargs="*"
)
argparser.add_argument("--datatype", help="Datatype", type=str, required=True)
argparser.add_argument("--timestamp", help="Timestamp", type=str, required=True)
argparser.add_argument("--configs", help="config file", type=str)
argparser.add_argument("--chan_maps", help="chan map", type=str)
argparser.add_argument("--metadata", help="metadata", type=str)
argparser.add_argument("--log", help="log file", type=str)
args = argparser.parse_args()

configs = TextDB(args.configs, lazy=True)
config_dict = configs.on(args.timestamp, system=args.datatype)["snakemake_rules"]["tier_raw"]

if "logging" in config_dict["options"]:
    log_config = config_dict["options"]["logging"]
    log_config = Props.read_from(log_config)
    if args.log is not None:
        Path(args.log).parent.mkdir(parents=True, exist_ok=True)
        log_config["handlers"]["file"]["filename"] = args.log
    logging.config.dictConfig(log_config)
    log = logging.getLogger(config_dict["options"].get("logger", "prod"))
else:
    if args.log is not None:
        Path(args.log).parent.makedir(parents=True, exist_ok=True)
        logging.basicConfig(level=logging.INFO, filename=args.log, filemode="w")
    log = logging.getLogger(__name__)

channel_dict = config_dict["inputs"]
hdf_settings = Props.read_from(config_dict["settings"])["hdf5_settings"]
blinding_settings = Props.read_from(config_dict["config"])

centroid = blinding_settings["centroid_in_keV"]  # keV
width = blinding_settings["width_in_keV"]  # keV

# list of all channels and objects in the raw file
all_channels = lh5.ls(args.input)

# list of Ge channels and SiPM channels with associated metadata
legendmetadata = LegendMetadata(args.metadata, lazy=True)
ged_channels = (
    legendmetadata.channelmap(args.timestamp).map("system", unique=False)["geds"].map("daq.rawid")
)
spms_channels = (
    legendmetadata.channelmap(args.timestamp).map("system", unique=False)["spms"].map("daq.rawid")
)
auxs_channels = (
    legendmetadata.channelmap(args.timestamp).map("system", unique=False)["auxs"].map("daq.rawid")
)
blsn_channels = (
    legendmetadata.channelmap(args.timestamp).map("system", unique=False)["bsln"].map("daq.rawid")
)
puls_channels = (
    legendmetadata.channelmap(args.timestamp).map("system", unique=False)["puls"].map("daq.rawid")
)

store = lh5.LH5Store()

# rows that need blinding
toblind = np.array([])

# first, loop through the Ge detector channels, calibrate them and look for events that should be blinded
for chnum in list(ged_channels):
    # skip Ge detectors that are anti-coincidence only or not able to be blinded for some other reason
    if ged_channels[chnum]["analysis"]["is_blinded"] is False:
        continue

    # load in just the daqenergy for now
    daqenergy, _ = store.read(f"ch{chnum}/raw/daqenergy", args.input)

    # read in calibration curve for this channel
    blind_curve = Props.read_from(args.blind_curve)[f"ch{chnum}"]["pars"]["operations"]

    # calibrate daq energy using pre existing curve
    daqenergy_cal = ne.evaluate(
        blind_curve["daqenergy_cal"]["expression"],
        local_dict=dict(daqenergy=daqenergy, **blind_curve["daqenergy_cal"]["parameters"]),
    )

    # figure out which event indices should be blinded
    toblind = np.append(
        toblind, np.nonzero(np.abs(np.asarray(daqenergy_cal) - centroid) <= width)[0]
    )

# remove duplicates
toblind = np.unique(toblind)

# total number of events (from last Ge channel loaded, should be same for all Ge channels)
allind = np.arange(len(daqenergy))

# gets events that should not be blinded
tokeep = allind[np.logical_not(np.isin(allind, toblind))]

# make some temp file to write the output to before renaming it
rng = np.random.default_rng()
rand_num = f"{rng.integers(0,99999):05d}"
temp_output = f"{args.output}.{rand_num}"
Path(temp_output).parent.mkdir(parents=True, exist_ok=True)

for channel in all_channels:
    try:
        chnum = int(channel[2::])
    except ValueError:
        # if this isn't an interesting channel, just copy it to the output file
        chobj, _ = store.read(channel, args.input, decompress=False)
        store.write_object(
            chobj,
            channel,
            lh5_file=temp_output,
            wo_mode="w",
            **hdf_settings,
        )
        continue

    if (
        (chnum not in list(ged_channels))
        and (chnum not in list(spms_channels))
        and (chnum not in list(auxs_channels))
        and (chnum not in list(blsn_channels))
        and (chnum not in list(puls_channels))
    ):
        # if this is a PMT or not included for some reason, just copy it to the output file
        chobj, _ = store.read(channel + "/raw", args.input, decompress=False)
        store.write_object(
            chobj,
            group=channel,
            name="raw",
            lh5_file=temp_output,
            wo_mode="w",
            **hdf_settings,
        )
        continue

    # the rest should be the Ge and SiPM channels that need to be blinded

    # read in all of the data but only for the unblinded events
    blinded_chobj, _ = store.read(channel + "/raw", args.input, idx=tokeep, decompress=False)

    # now write the blinded data for this channel
    store.write_object(
        blinded_chobj,
        group=channel,
        name="raw",
        lh5_file=temp_output,
        wo_mode="w",
        **hdf_settings,
    )

# rename the temp file
Path(args.output).parent.mkdir(parents=True, exist_ok=True)
Path(temp_output).rename(args.output)
