import pytest
from pathlib import Path

import pygama.lgdo as lgdo
import numpy as np
import json
import os

from pygama.data_trimmer.data_trimmer import data_trimmer


config_dir = Path(__file__).parent / "test_data_trimmer_configs"

# check that packet indexes match in verification test 
def test_data_trimmer_packet_ids(lgnd_test_data):

    # Set up I/O files, including config 
    raw_file = lgnd_test_data.get_path("lh5/prod-ref-l200/generated/tier/raw/cal/p01/r014/l60-p01-r014-cal-20220716T104550Z-tier_raw.lh5")
    dsp_config = f"{config_dir}/data_trimmer_config.json"

    presummed_file = raw_file.replace("l60-p01-r014-cal-20220716T104550Z-tier_raw.lh5", "l60-p01-r014-cal-20220716T104550Z-tier_raw_presummed.lh5")
    windowed_file = raw_file.replace("l60-p01-r014-cal-20220716T104550Z-tier_raw.lh5", "l60-p01-r014-cal-20220716T104550Z-tier_raw_windowed.lh5")

    data_trimmer(raw_file, windowed_file, presummed_file, dsp_config)

    lh5_tables = lgdo.ls(raw_file)
    # check if group points to raw data; sometimes 'raw' is nested, e.g g024/raw
    for i, tb in enumerate(lh5_tables):
        if "raw" not in tb and lgdo.ls(raw_file, f"{tb}/raw"):
            lh5_tables[i] = f"{tb}/raw"
        elif not lgdo.ls(raw_file, tb):
            del lh5_tables[i]

    sto = lgdo.LH5Store()

    for raw_group in lh5_tables: 
        raw_packet_ids = sto.read_object(str(raw_group) + "/packet_id", raw_file)
        presummed_packet_ids = sto.read_object(str(raw_group) + "/packet_id", presummed_file)
        windowed_packet_ids = sto.read_object(str(raw_group) + "/packet_id", windowed_file)

        assert np.array_equal(raw_packet_ids[0].nda, presummed_packet_ids[0].nda)
        assert np.array_equal(raw_packet_ids[0].nda, windowed_packet_ids[0].nda)


# check that packet indexes match in verification test 
def test_data_trimmer_waveform_lengths(lgnd_test_data):

    # Set up I/O files, including config 
    raw_file = lgnd_test_data.get_path("lh5/prod-ref-l200/generated/tier/raw/cal/p01/r014/l60-p01-r014-cal-20220716T105236Z-tier_raw.lh5")
    presummed_file = raw_file.replace("l60-p01-r014-cal-20220716T105236Z-tier_raw.lh5", "l60-p01-r014-cal-20220716T105236Z-tier_raw_presummed.lh5")
    windowed_file = raw_file.replace("l60-p01-r014-cal-20220716T105236Z-tier_raw.lh5", "l60-p01-r014-cal-20220716T105236Z-tier_raw_windowed.lh5")

    dsp_config = '''
    {
        "outputs" : [ "presummed" ],
        "processors" : {            
            "presummed": {
                "function": "presum",
                "module": "pygama.dsp.processors",
                "args": ["waveform", "presummed(len(waveform)/4, 'f')"],
                "unit": "ADC"
            }
        }
    }
    '''

    data_trimmer(raw_file, windowed_file, presummed_file, dsp_config)

    lh5_tables = lgdo.ls(raw_file)
    # check if group points to raw data; sometimes 'raw' is nested, e.g g024/raw
    for i, tb in enumerate(lh5_tables):
        if "raw" not in tb and lgdo.ls(raw_file, f"{tb}/raw"):
            lh5_tables[i] = f"{tb}/raw"
        elif not lgdo.ls(raw_file, tb):
            del lh5_tables[i]

    if isinstance(dsp_config, str) and dsp_config.endswith(".json"):
        f = open(dsp_config)
        jsonfile = json.load(f)
        f.close()
    # If we get a string that is in the correct format as a json file
    elif isinstance(dsp_config, str):
        jsonfile = json.loads(dsp_config)
    # Or we could get a dict as the config
    elif isinstance(dsp_config, dict):
        jsonfile = dsp_config

    # Read in the presummed rate from the config file to modify the clock rate later 
    presum_rate_string = jsonfile["processors"]["presummed"]["args"][1]
    presum_rate_start_idx = presum_rate_string.find("/")+1
    presum_rate_end_idx = presum_rate_string.find(",")
    presum_rate = int(presum_rate_string[presum_rate_start_idx:presum_rate_end_idx])

    # This needs to be overwritten with the correct windowing values set in data_trimmer.py 
    window_start_index = 1000
    window_end_index = 1000


    sto = lgdo.LH5Store()

    for raw_group in lh5_tables:

        raw_packet_waveform_values = sto.read_object(str(raw_group) + "/waveform/values", raw_file)
        presummed_packet_waveform_values = sto.read_object(str(raw_group) +"/waveform/values", presummed_file)
        windowed_packet_waveform_values = sto.read_object(str(raw_group) +"/waveform/values", windowed_file)

        assert len(raw_packet_waveform_values[0].nda[0]) == presum_rate*len(presummed_packet_waveform_values[0].nda[0])
        assert len(raw_packet_waveform_values[0].nda[0]) == len(windowed_packet_waveform_values[0].nda[0])+window_start_index+window_end_index

def test_data_trimmer_file_size_decrease(lgnd_test_data):
    # Set up I/O files, including config 
    raw_file = lgnd_test_data.get_path("lh5/LDQTA_r117_20200110T105115Z_cal_geds_raw.lh5")
    dsp_config = f"{config_dir}/data_trimmer_config.json"

    presummed_file = raw_file.replace("LDQTA_r117_20200110T105115Z_cal_geds_raw.lh5", "LDQTA_r117_20200110T105115Z_cal_geds_raw_presummed.lh5")
    windowed_file = raw_file.replace("LDQTA_r117_20200110T105115Z_cal_geds_raw.lh5", "LDQTA_r117_20200110T105115Z_cal_geds_raw_windowed.lh5")

    data_trimmer(raw_file, windowed_file, presummed_file, dsp_config)

    assert os.path.getsize(presummed_file) < os.path.getsize(raw_file)
    assert os.path.getsize(windowed_file) < os.path.getsize(raw_file)

