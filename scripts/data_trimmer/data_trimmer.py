from __future__ import annotations

import pygama.lgdo as lgdo
import numpy as np
from pygama.dsp.processing_chain import build_processing_chain as bpc
import json

def data_trimmer(raw_file: str, windowed_file: str, presummed_file: str, dsp_config: str | dict) -> None:
    """
    Takes in a raw file, and returns two files, one containing the presummed waveform, and the other 
    containing a windowed waveform.

    Parameters 
    ----------
    raw_file
        A raw lh5 file to window and presum
    windowed_file
        Name of an lh5 file that will contain windowed waveforms 
    presummed_file
        Name of an lh5 file that will contain presummed waveforms 
    dsp_config 
        Either the path to the dsp config file to use, or a dictionary of config
    
    Notes 
    ----- 
    The windowing indices need to be set inside this file, set window_start_index and window_end_index appropriately
    """
    # In the future, put a build_processing_chain-esque dsp and write statement into the read_chunk() loop part of build_raw

    # Read in the raw file
    sto = lgdo.LH5Store()
    lh5_tables = lgdo.ls(raw_file)

    # check if group points to raw data; sometimes 'raw' is nested, e.g g024/raw
    for i, tb in enumerate(lh5_tables):
        if "raw" not in tb and lgdo.ls(raw_file, f"{tb}/raw"):
            lh5_tables[i] = f"{tb}/raw"
        elif not lgdo.ls(raw_file, tb):
            del lh5_tables[i]

    if len(lh5_tables) == 0:
        raise RuntimeError(f"could not find any valid LH5 table in {f_raw}")

    # Grab the constants from the dsp config 
    if isinstance(dsp_config, str) and dsp_config.endswith(".json"):
        f = open(dsp_config)
        dsp_dict = json.load(f)
        f.close()
    # If we get a string that is in the correct format as a json file
    elif isinstance(dsp_config, str):
        dsp_dict = json.loads(dsp_config)
    # Or we could get a dict as the config
    elif isinstance(dsp_config, dict):
        dsp_dict = dsp_config


    # Read in the presummed rate from the config file to modify the clock rate later 
    presum_rate_string = dsp_dict["processors"]["presummed"]["args"][1]
    presum_rate_start_idx = presum_rate_string.find("/")+1
    presum_rate_end_idx = presum_rate_string.find(",")
    presum_rate = int(presum_rate_string[presum_rate_start_idx:presum_rate_end_idx])

    # Need to override this value to get the correct waveform window
    window_start_index = int(1000)
    window_end_index = int(1000)

    # loop over the tables in the file, allows for the case of multiple channels per file 
    for raw_group in lh5_tables:        
        # Read in the table that we will modify
        raw_table, _ = sto.read_object(raw_group, raw_file)
       
        # Execute the processing chain
        proc_chain, mask, dsp_out = bpc(raw_table, dsp_dict)
        proc_chain.execute()

        # Write the windowed waveform to a file 
        # Overwrite the waveform values with a windowed view of the same raw_table waveform values 
        # This is more memory efficient than using build_processing_chain to allocate a new array
        raw_table["waveform"]["values"].nda =  raw_table["waveform"]["values"].nda[ : , window_start_index:-window_end_index]

        # For the windowed waveform, change the t0 by the windowing (t0+start_index*dt?) or units of (t0/dt+start_index) in samples
        raw_table["waveform"]["t0"].nda /= raw_table["waveform"]["dt"].nda * float(presum_rate) # undo the clock rate transform from the presummed wf
        raw_table["waveform"]["t0"].nda += float(window_start_index)
        raw_table["waveform"]["t0"].attrs['units'] = "samples" # Changed the units to samples

        sto.write_object(raw_table, "raw", windowed_file, group=raw_group.split("/raw")[0], wo_mode="a")

        # Write the presummed waveform to file 
        # Overwrite the waveform values with those from the processed waveform
        raw_table["waveform"]["values"].nda = dsp_out["presummed"].nda

        # For the presummed waveform, overwrite the dt:
        raw_table["waveform"]["dt"].nda *= float(presum_rate)
        
        # Write to raw files
        sto.write_object(raw_table, "raw", presummed_file, group=raw_group.split("/raw")[0], wo_mode="a")