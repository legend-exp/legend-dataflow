from .utils import *
from .patterns import *
import json

def gen_file_db_config(setup):
    db_config = {
        "data_dir" : "",
        "tier_dirs":{
            "raw": "/raw",
            "dsp": "/dsp",
            "hit": "/hit",
            "tcm": "/tcm",
            "evt": "/evt"
        },
        "file_format": {
            "raw": get_pattern_tier(setup,"raw"),
            "dsp": get_pattern_tier(setup,"dsp"),    # need to remove data_dir from path?  
            "hit": get_pattern_tier(setup,"hit"),    # can this handle different dirs as raw data will not be in same path in future
            "evt": get_pattern_tier(setup,"evt"),
            "tcm": get_pattern_tier(setup,"tcm")
            },
        "table_format":{
            "raw": "ch{ch:03d}/raw",                # this should be checked from raw config and then applied to higher tiers
            "dsp": "ch{ch:03d}/dsp",
            "hit": "ch{ch:03d}/hit",
            "evt": "{grp}/evt",
            "tcm": "hardware_tcm"
        }
    }

setup = snakemake.params.setup

with open(output[0], "w") as w:
    json.dump(gen_file_db_config(setup) ,w ,indent=4)