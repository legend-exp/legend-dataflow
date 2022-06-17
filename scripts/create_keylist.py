#import snakemake as smk
#import os, re, glob

#from utils import *
from util.FileKey import *
import glob

keypart = snakemake.wildcards.keypart
setup = snakemake.params.setup

"""
d = parse_keypart(keypart)


tier0_pattern = get_pattern_evts(setup, "daq")
tier0_pattern_rx = re.compile(smk.io.regex(tier0_pattern))
fn_glob_pattern = smk.io.expand(tier0_pattern, experiment = d["experiment"], period = d["period"], run = d["run"], datatype = d["datatype"], timestamp = d["timestamp"])[0]

files = glob.glob(fn_glob_pattern)

keys = []
for f in files:
    m  = tier0_pattern_rx.match(f)
    if m is not None:

        d = m.groupdict()
        key = smk.io.expand(f"{key_pattern()}", experiment = d["experiment"], period = d["period"], run = d["run"], datatype  = d["datatype"], timestamp = d["timestamp"])[0]
        keys.append(key)

with open(snakemake.output[0], 'w') as f:
    for key in keys:
        f.write(f"{key}\n")
"""

key = FileKey.parse_keypart('-l200')
fn_glob_pattern = key.get_path_from_filekey(get_pattern_evts(setup, "daq"))[0]
files = glob.glob(fn_glob_pattern)
keys = []
for f in files:
    key = FileKey.get_filekey_from_pattern(f,get_pattern_evts(setup, "daq"))
    keys.append(key.name)

with open(snakemake.output[0], 'w') as f:
    for key in keys:
        f.write(f"{key}\n")