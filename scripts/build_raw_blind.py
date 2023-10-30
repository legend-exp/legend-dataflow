"""
This script takes in raw data, applies the calibration to the daqenergy
and uses this to blind the data in a window of Qbb +- 25 keV. It copies over all
channels in a raw file, removing those events that fall within the ROI for Ge detectors
that have a daqenergy calibration curve and are not anti-coincidence only (AC).

In the Snakemake dataflow, this script only runs if the checkfile is found on disk,
but this is controlled by the Snakemake flow (presumably an error is thrown if the file
is not found). This script itself does not check for the existence of such a file.
"""

import argparse
import logging
import os
import pathlib
import numexpr as ne
import numpy as np
from legendmeta import LegendMetadata
from legendmeta.catalog import Props
import lgdo.lh5_store as lh5

argparser = argparse.ArgumentParser()
argparser.add_argument("input", help="input file", type=str)
argparser.add_argument("output", help="output file", type=str)
argparser.add_argument("--blind_curve", help="blinding curves file", type=str, required=True)
argparser.add_argument("--datatype", help="Datatype", type=str, required=True)
argparser.add_argument("--timestamp", help="Timestamp", type=str, required=True)
argparser.add_argument("--configs", help="config file", type=str)
argparser.add_argument("--chan_maps", help="chan map", type=str)
argparser.add_argument("--log", help="log file", type=str)
args = argparser.parse_args()

os.makedirs(os.path.dirname(args.log), exist_ok=True)
logging.basicConfig(level=logging.INFO, filename=args.log, filemode="w")

pathlib.Path(os.path.dirname(args.output)).mkdir(parents=True, exist_ok=True)

Qbb = 2039.061 # keV
ROI = 25.0 # keV

# list of all channels and objects in the raw file
all_channels = lh5.ls(args.input)

# list of just germanium channels with associated metadata
# I'm not sure if this is supposed to be "daq.fcid" or "daq.fc_channel" 
# (from a recent pull that renamed it) or "daq.rawid", which is what I've chosen for now
chmap = LegendMetadata(path=args.chan_maps)
ged_channels = chmap.channelmap(args.timestamp).map("system", unique=False)["geds"].map("daq.rawid")

store = lh5.LH5Store()

for channel in all_channels:
    try:
        chnum = int(channel[2::])
    except ValueError:
        # if this isn't an interesting channel, just copy it to the output file
        chobj, _ = store.read_object(channel, args.input)
        store.write_object(chobj, channel, args.output, wo_mode='overwrite')
        continue
    
    if chnum not in list(ged_channels):
        # if this is a SiPM or Ge not included for some reason, just copy it to the output file
        chobj, _ = store.read_object(channel+'/raw', args.input)
        store.write_object(chobj, group=channel, name='raw', lh5_file=args.output, wo_mode='overwrite')
        continue

    if ged_channels[chnum]["analysis"]["usability"] == 'ac':
        # if this Ge is to be used for anti-coincidence only, it will not have a blinding calibration
        # (or at least it should not be blinded) so just copy it to the output file
        chobj, _ = store.read_object(channel+'/raw', args.input)
        store.write_object(chobj, group=channel, name='raw', lh5_file=args.output, wo_mode='overwrite')
        continue

    # the rest should be the Ge channels that need to be blinded

    # load in just the daqenergy for now
    daqenergy, _  = store.read_object(channel+'/raw/daqenergy', args.input)

    # read in calibration curve for this channel
    blind_curve = Props.read_from(args.blind_curve)[channel]

    # calibrate daq energy using pre existing curve
    daqenergy_cal = ne.evaluate(
        blind_curve["daqenergy_cal"]["expression"],
        local_dict=dict(daqenergy=daqenergy, **blind_curve["daqenergy_cal"]["parameters"]),
    )

    # figure out which events should be kept unblinded
    tokeep = np.nonzero(np.abs(np.asarray(daqenergy_cal) - Qbb) > ROI)[0]

    # read in all of the data but only for the unblinded events
    # this says 'idx empty after culling.' but it seems to work?
    # I made a pull request for lgdo to fix this bug.
    blinded_chobj, _  = store.read_object(channel+'/raw', args.input, idx=tokeep)

    # now write the blinded data for this channel
    store.write_object(blinded_chobj, group=channel, name='raw', lh5_file=args.output, wo_mode='overwrite')
    
# think this was probably for testing
# pathlib.Path(args.output).touch()
