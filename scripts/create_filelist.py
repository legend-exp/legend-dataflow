import snakemake as smk
import os, re, glob, json

from util.FileKey import *

setup = snakemake.params.setup

file_selection = snakemake.wildcards.label[:3]
tier = snakemake.wildcards.tier
keypart=snakemake.wildcards.label[3:]#.replace("all","")
search_pattern = snakemake.params.search_pattern

ignore_keys = []
if snakemake.params.configs:
    configs=snakemake.params.configs
    ignored_keyslist = os.path.join(configs, 'ignore_keys.keylist')
    if os.path.isfile(ignored_keyslist):
        with open(ignored_keyslist, 'r') as f:
            ignore_keys = f.read().splitlines()
    else:
        print("no ignore_keys.keylist file found")

    analysis_runs_file = os.path.join(configs, 'analysis_runs.json')
    if os.path.isfile(analysis_runs_file):
        with open(analysis_runs_file, 'r') as f:
            analysis_runs = json.load(f)
    else:
        print("no analysis_runs file found")

key = FileKey.parse_keypart(keypart)

item_list = []
for item in key:
    if "_" in item:
        item = item.split("_")
    if isinstance(item, list):
        item_list.append(item)
    else:
        item_list.append([item])

filekeys=[]
for i in item_list[0]:
    for j in item_list[1]:
        for k in item_list[2]:
            for l in item_list[3]:
                for m in item_list[4]:
                    filekeys.append(FileKey(i,j,k,l,m))

keys = []
for key in filekeys:
    fn_glob_pattern = key.get_path_from_filekey(search_pattern)[0]
    files = glob.glob(fn_glob_pattern)
    for f in files:
        key = FileKey.get_filekey_from_pattern(f,search_pattern)
        if key.name in ignore_keys:
            pass
        else:
            if file_selection == 'all':
                keys.append(key.name)
            elif file_selection == "sel":
                if analysis_runs == "all":
                    keys.append(key.name)
                else:
                    if key.period in analysis_runs:
                        if key.run in analysis_runs[key.period] or analysis_runs[key.period]=="all":
                            keys.append(key.name)
            else:
                raise ValueError("unknown file selection")
keys = sorted(keys)


filenames = FileKey.tier_files(setup, keys, tier)
with open(snakemake.output[0], 'w') as f:
    for fn in filenames:
        f.write(f"{fn}\n")

