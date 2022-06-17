import snakemake as smk
import os, re, glob
import json
from .utils import *
import argparse
from .FileKey import *
from .patterns import *

class pars_key_resolve():
    
    name_dict = {"cal":["dsp","hit"], 'lar':['dsp','hit']}
    
    def __init__(self, valid_from, category, apply):
        self.valid_from = valid_from
        self.category = category
        self.apply = apply
    
    def get_json(self):
        return json.dumps(self.__dict__)
    
    @classmethod
    def from_filekey(cls, filekey):
        return cls(filekey.timestamp, "all", filekey.get_file(f'{par_pattern()}.json', pars_key_resolve.name_dict))
    
    @staticmethod
    def write_to_jsonl(file_names, path):
        with open(path, "w") as of:
            for file_name in file_names:
                of.write(f'{file_name.get_json()}\n')
                
    @staticmethod
    def match_keys(key1, key2):
        if key1.experiment == key2.experiment and key1.period == key2.period and key1.run == key2.run and key1.datatype == key2.datatype:
            if key1.get_unix_timestamp() < key2.get_unix_timestamp():
                return key1
            else:
                return key2
        else:
            return key2
    
    @staticmethod
    def generate_par_keylist(keys):
        keylist = []
        keylist.append(keys[0])
        for key in keys[1:]:
            matched_key = pars_key_resolve.match_keys(keylist[-1], key)
            if matched_key not in keylist:
                keylist.append(matched_key)
            else:
                pass
        return keylist
    
    @staticmethod    
    def match_entries(entry1, entry2):
        datatype2 = FileKey.get_filekey_from_filename(entry2.apply[0]).datatype
        for entry in entry1.apply:
            if FileKey.get_filekey_from_filename(entry).datatype == datatype2:
                pass
            else:
                entry2.apply.append(entry)
    
    @staticmethod
    def match_all_entries(entrylist):
        out_list = []
        out_list.append(pars_key_resolve.from_filekey(entrylist[0]))
        for entry in entrylist[1:]:
            new_entry = pars_key_resolve.from_filekey(entry)
            pars_key_resolve.match_entries(out_list[-1],new_entry)
            out_list.append(new_entry)
        return out_list
    
    @staticmethod
    def get_keys(keypart, setup):
        d = FileKey.parse_keypart(keypart)
        tier0_pattern = get_pattern_evts(setup, "daq")
        tier0_pattern_rx = re.compile(smk.io.regex(tier0_pattern))
        fn_glob_pattern = smk.io.expand(tier0_pattern, **d._asdict())[0]
        files = glob.glob(fn_glob_pattern)
        keys = []
        for f in files:
            m  = tier0_pattern_rx.match(f)
            if m is not None:
                d = m.groupdict()
                key = FileKey(**d)
                keys.append(key)
        return keys
    
    @staticmethod
    def write_par_catalog(setup, keypart, filename):
        if isinstance(keypart, str):
            keypart = [keypart]
        keylist = []
        for keypar in keypart:
            keylist += pars_key_resolve.get_keys(keypar, setup)
        keys = sorted(keylist, key=FileKey.get_unix_timestamp)
        keylist = pars_key_resolve.generate_par_keylist(keys)
        entrylist = pars_key_resolve.match_all_entries(keylist)
        pars_key_resolve.write_to_jsonl(entrylist, filename)
        