"""
This module creates the validity files used for determining the time validity of data
"""

import json
import re
import warnings
from pathlib import Path

import snakemake as smk
import yaml

from .FileKey import FileKey, ProcessingFileKey
from .patterns import par_validity_pattern


class ParsKeyResolve:

    def __init__(self, valid_from, category, apply):
        self.valid_from = valid_from
        self.category = category
        self.mode = "reset"
        self.apply = apply

    def __str__(self):
        return f"{self.__dict__}"

    def get_json(self):
        return json.dumps(self.__dict__)

    @classmethod
    def from_filekey(cls, filekey, name_dict):
        return cls(
            filekey.timestamp,
            "all",
            filekey.get_path_from_filekey(
                par_validity_pattern(), processing_step=name_dict, ext="yaml"
            ),
        )

    @staticmethod
    def write_to_jsonl(file_names, path):
        with Path(path).open("w") as of:
            for file_name in file_names:
                of.write(f"{file_name.get_json()}\n")

    @staticmethod
    def write_to_yaml(file_names, path):
        with Path(path).open("w") as of:
            yaml.dump([file_name.__dict__ for file_name in file_names], of, sort_keys=False)

    @staticmethod
    def match_keys(key1, key2):
        if (
            key1.experiment == key2.experiment
            and key1.period == key2.period
            and key1.run == key2.run
            and key1.datatype == key2.datatype
        ):
            if key1.get_unix_timestamp() < key2.get_unix_timestamp():
                return key1
            else:
                return key2
        else:
            return key2

    @staticmethod
    def generate_par_keylist(keys):
        keylist = []
        keys = sorted(keys, key=FileKey.get_unix_timestamp)
        keylist.append(keys[0])
        for key in keys[1:]:
            matched_key = ParsKeyResolve.match_keys(keylist[-1], key)
            if matched_key not in keylist:
                keylist.append(matched_key)
            else:
                pass
        return keylist

    @staticmethod
    def match_entries(entry1, entry2):
        datatype2 = ProcessingFileKey.get_filekey_from_filename(entry2.apply[0]).datatype
        for entry in entry1.apply:
            if ProcessingFileKey.get_filekey_from_filename(entry).datatype == datatype2:
                pass
            else:
                entry2.apply.append(entry)

    @staticmethod
    def match_all_entries(entrylist, name_dict):
        out_list = []
        out_list.append(ParsKeyResolve.from_filekey(entrylist[0], name_dict))
        for entry in entrylist[1:]:
            new_entry = ParsKeyResolve.from_filekey(entry, name_dict)
            ParsKeyResolve.match_entries(out_list[-1], new_entry)
            out_list.append(new_entry)
        return out_list

    @staticmethod
    def get_keys(keypart, search_pattern):
        d = FileKey.parse_keypart(keypart)
        if Path(search_pattern).suffix == ".*":
            search_pattern = Path(search_pattern).with_suffix(".{ext}")
            wildcard_dict = dict(ext="*", **d._asdict())
        else:
            wildcard_dict = d._asdict()
        try:
            tier_pattern_rx = re.compile(smk.io.regex_from_filepattern(str(search_pattern)))
        except AttributeError:
            tier_pattern_rx = re.compile(smk.io.regex(str(search_pattern)))
        fn_glob_pattern = smk.io.expand(search_pattern, **wildcard_dict)[0]
        p = Path(fn_glob_pattern)
        parts = p.parts[p.is_absolute() :]
        files = Path(p.root).glob(str(Path(*parts)))
        keys = []
        for f in files:
            m = tier_pattern_rx.match(str(f))
            if m is not None:
                d = m.groupdict()
                if "ext" in d:
                    d.pop("ext")
                key = FileKey(**d)
                keys.append(key)
        return keys

    @staticmethod
    def get_par_catalog(keypart, search_patterns, name_dict):
        if isinstance(keypart, str):
            keypart = [keypart]
        if isinstance(search_patterns, (str, Path)):
            search_patterns = [search_patterns]
        keylist = []
        for search_pattern in search_patterns:
            for keypar in keypart:
                keylist += ParsKeyResolve.get_keys(keypar, search_pattern)
        if len(keylist) != 0:
            keys = sorted(keylist, key=FileKey.get_unix_timestamp)
            keylist = ParsKeyResolve.generate_par_keylist(keys)

            entrylist = ParsKeyResolve.match_all_entries(keylist, name_dict)
        else:
            msg = "No Keys found"
            warnings.warn(msg, stacklevel=0)
            entrylist = [ParsKeyResolve("00000000T000000Z", "all", [])]
        return entrylist
