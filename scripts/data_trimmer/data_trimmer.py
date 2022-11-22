from __future__ import annotations

import pygama.lgdo as lgdo
import numpy as np
from pygama.dsp.processing_chain import build_processing_chain as bpc
import json
import copy


def data_trimmer(raw_file: str, trimmed_file: str, dsp_config: str | dict) -> None:
    """
    Takes in a raw file, and returns two files, one containing the presummed waveform, and the other 
    containing a windowed waveform.

    Parameters 
    ----------
    raw_file
        A raw lh5 file to window and presum
    trimmed_file
        Name of an lh5 file that will contain windowed and trimmed waveforms, under the group ``.../raw/presummed_waveform`` and ``.../raw/windowed_waveform``
    dsp_config 
        Either the path to the dsp config file to use, or a dictionary of config
    
    Notes 
    ----- 
    The windowing indices need to be set inside this file, set window_start_index and window_end_index appropriately
    The presummed waveform is unnormalized! The normalization is stored inside a table of the same length as t0 in waveform as "presum_len" 
    """
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
        # Just manipulate the existing raw table so we don't have to allocate new memory 
        raw_table["windowed_waveform"] = raw_table.pop('waveform')
        raw_table.update_datatype()

        # Overwrite the waveform values with a windowed view of the same raw_table waveform values 
        # This is more memory efficient than using build_processing_chain to allocate a new array
        raw_table["windowed_waveform"]["values"].nda =  raw_table["windowed_waveform"]["values"].nda[ : , window_start_index:-window_end_index]
        
        original_t0_units = raw_table["windowed_waveform"]["t0"].attrs['units']

        # For the windowed waveform, change the t0 by the windowing to units of (t0/dt+start_index) in samples
        raw_table["windowed_waveform"]["t0"].nda /= raw_table["windowed_waveform"]["dt"].nda 
        raw_table["windowed_waveform"]["t0"].nda += float(window_start_index)
        raw_table["windowed_waveform"]["t0"].attrs['units'] = "samples" # Changed the units to samples

        # Create a new table that contains the presummed waveform WaveformTable
        presummed_table = lgdo.WaveformTable(size= len(raw_table["windowed_waveform"]["t0"].nda), t0 = 0, t0_units=copy.deepcopy(raw_table["windowed_waveform"]["t0"].attrs['units']),  dt = 0, dt_units= copy.deepcopy(raw_table["windowed_waveform"]["dt"].attrs['units']), wf_len= int((len(raw_table["windowed_waveform"]["values"].nda[0])+window_start_index+window_end_index)/presum_rate) , dtype= np.uint32 )
        raw_table.add_field("presummed_waveform", presummed_table, use_obj_size = True)
        
        # Write the presummed waveform to file 
        # Overwrite the waveform values with those from the processed waveform
        raw_table["presummed_waveform"]["values"].nda = dsp_out["presummed"].nda

        # Change the t0 back to the original units and values
        raw_table["presummed_waveform"]["t0"].nda -= float(window_start_index)
        raw_table["presummed_waveform"]["t0"].nda *= raw_table["presummed_waveform"]["dt"].nda 
        raw_table["presummed_waveform"]["t0"].attrs['units'] = str(original_t0_units) 
        
        # For the presummed waveform, overwrite the dt:
        raw_table["presummed_waveform"]["dt"].nda = float(presum_rate)*raw_table["windowed_waveform"]["dt"].nda
        
        # Write to raw files
        sto.write_object(raw_table, "raw", trimmed_file, group=raw_group.split("/raw")[0], wo_mode="a")