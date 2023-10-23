# ruff: noqa: F821, T201

import glob
import os

from util.FileKey import *

ignore_keys = []
if snakemake.input:
    configs = snakemake.input[0]
    ignored_keyslist = os.path.join(configs, "ignore_keys.keylist")
    if os.path.isfile(ignored_keyslist):
        with open(ignored_keyslist) as f:
            ignore_keys = f.read().splitlines()
    else:
        print("no ignore_keys.keylist file found")


keypart = snakemake.wildcards.keypart
setup = snakemake.params.setup
search_pattern = snakemake.params.search_pattern

key = FileKey.parse_keypart(keypart)

item_list = []
for item in key:
    if "_" in item:
        _item = item.split("_")
    if isinstance(_item, list):
        item_list.append(_item)
    else:
        item_list.append([_item])

filekeys = []
for i in item_list[0]:
    for j in item_list[1]:
        for k in item_list[2]:
            for i2 in item_list[3]:
                for j2 in item_list[4]:
                    filekeys.append(FileKey(i, j, k, i2, j2))

keys = []
for key in filekeys:
    fn_glob_pattern = key.get_path_from_filekey(search_pattern)[0]
    files = glob.glob(fn_glob_pattern)

    for f in files:
        _key = FileKey.get_filekey_from_pattern(f, search_pattern)
        if _key.name in ignore_keys:
            pass
        else:
            keys.append(_key.name)

with open(snakemake.output[0], "w") as f:
    for key in keys:
        f.write(f"{key}\n")
