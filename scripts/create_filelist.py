import snakemake as smk
import os, re, glob

from utils import *

setup = snakemake.params.setup

dataset_file = snakemake.input[0]
tier = snakemake.wildcards.tier
keypart = snakemake.wildcards.label

if "channels" in keypart:
    name = f'{tier}_pars'
    
    d = parse_channel_keypart(keypart)

    par_pattern = get_pattern_pars_tmp_channel(setup,tier,name)
    par_pattern_rx = re.compile(smk.io.regex(par_pattern))
    fn_glob_pattern = smk.io.expand(par_pattern, experiment = d["experiment"], 
                                    period = d["period"], run = d["run"], 
                                    datatype = d["datatype"], channel=d["channel"])[0]

    filenames = []
    with open(dataset_file) as f:
        for channel in f:
            file = fn_glob_pattern.replace('channels', channel)
            filenames.append(file)

    with open(snakemake.output[0], 'w') as f:
        for fn in filenames:
            f.write(f"{fn}\n")
        #f.write("/lfs/l1/legend/legend-prodenv/prod-usr/ggmarsh-l200_prototype-v01/generated/tmp/par/dsp/cal/p001/r0001/dsp_pars/1/dsp_pars-l200-p001-r0001-1-dsp.json")
else:
    filenames = tier_files(setup, dataset_file, tier)
    with open(snakemake.output[0], 'w') as f:
        for fn in filenames:
            f.write(f"{fn}\n")
        #f.write("/lfs/l1/legend/legend-prodenv/prod-usr/ggmarsh-l200_prototype-v01/generated/par/dsp/cal/p001/r0001/dsp_pars-l200-p001-r0001-cal-dsp.json")
