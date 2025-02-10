"""
This module creates the validity files used for determining the time validity of data
"""

import re
import warnings
from pathlib import Path

from dbetto import time
from dbetto.catalog import Catalog

from .FileKey import FileKey, ProcessingFileKey, regex_from_filepattern
from .patterns import par_validity_pattern


class ParsKeyResolve(Catalog):
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

    @classmethod
    def entry_from_filekey(cls, filekey, name_dict):
        return cls.Entry(
            time.unix_time(filekey.timestamp),
            filekey.get_path_from_filekey(
                par_validity_pattern(), processing_step=name_dict, ext="yaml"
            ),
        )

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
        datatype2 = ProcessingFileKey.get_filekey_from_filename(entry2.file[0]).datatype
        for entry in entry1.file:
            if ProcessingFileKey.get_filekey_from_filename(entry).datatype == datatype2:
                pass
            else:
                entry2.file.append(entry)

    @staticmethod
    def match_all_entries(entrylist, name_dict):
        out_list = []
        out_list.append(ParsKeyResolve.entry_from_filekey(entrylist[0], name_dict))
        for entry in entrylist[1:]:
            new_entry = ParsKeyResolve.entry_from_filekey(entry, name_dict)
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

        tier_pattern_rx = re.compile(regex_from_filepattern(str(search_pattern)))
        key = FileKey.get_filekey_from_pattern(search_pattern, search_pattern)
        fn_glob_pattern = key.get_path_from_filekey(search_pattern, **wildcard_dict)[0]
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

    @classmethod
    def get_par_catalog(cls, keypart, search_patterns, name_dict):
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
            fk = FileKey("l200", "p00", "r000", "cal", "20230101T000000Z")
            entrylist = [ParsKeyResolve.entry_from_filekey(fk, name_dict)]
        return Catalog({"all": entrylist})
