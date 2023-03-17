import snakemake as smk
import os, re, glob

from util.FileKey import *

setup = snakemake.params.setup

tier = snakemake.wildcards.tier
keypart=snakemake.wildcards.label.replace("all","")
search_pattern = snakemake.params.search_pattern

ignore_keys = []
if snakemake.input:
    configs=snakemake.params.configs
    ignored_keyslist = os.path.join(configs, 'ignore_keys.keylist')
    if os.path.isfile(ignored_keyslist):
        with open(ignored_keyslist, 'r') as f:
            ignore_keys = f.read().splitlines()
    else:
        print("no ignore_keys.keylist file found")

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
            keys.append(key.name)


filenames = FileKey.tier_files(setup, keys, tier)
with open(snakemake.output[0], 'w') as f:
    for fn in filenames:
        f.write(f"{fn}\n")

