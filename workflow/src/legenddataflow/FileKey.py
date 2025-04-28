"""
This module contains classes to convert between keys and files using the patterns defined in patterns.py
"""

import re
import string
from collections import namedtuple
from itertools import product
from pathlib import Path

from dbetto.time import unix_time

from .patterns import (
    full_channel_pattern_with_extension,
    get_pattern_tier,
    key_pattern,
    processing_pattern,
)

# key_pattern -> key
#


def regex_from_filepattern(filepattern):
    f = []
    wildcards = []
    last = 0
    for match in re.compile(r"\{(?P<name>[\w]+)\}").finditer(filepattern):
        f.append(re.escape(filepattern[last : match.start()]))
        wildcard = match.group("name")
        if wildcard in wildcards:
            f.append(f"(?P={wildcard})")
        else:
            wildcards.append(wildcard)
            if wildcard == "ext":
                f.append(
                    f"(?P<{wildcard}>.*)"
                )  # this means ext will capture everything after 1st dot
            else:
                f.append(f"(?P<{wildcard}>" + r"[^\.\/]+)")
        last = match.end()
    f.append(re.escape(filepattern[last:]))
    f.append("$")
    return "".join(f)


class FileKey(
    namedtuple("FileKey", ["experiment", "period", "run", "datatype", "timestamp"])
):
    __slots__ = ()

    re_pattern = "(-(?P<experiment>[^-]+)(\\-(?P<period>[^-]+)(\\-(?P<run>[^-]+)(\\-(?P<datatype>[^-]+)(\\-(?P<timestamp>[^-]+))?)?)?)?)?$"
    key_pattern = key_pattern()

    @property
    def name(self):
        return f"{self.experiment}-{self.period}-{self.run}-{self.datatype}-{self.timestamp}"

    @property
    def key(self):
        return f"-{self.name}"

    def _list(self):
        return list(self)

    @property
    def __str__(self):
        return self.name

    @classmethod
    def from_string(cls, key_string):
        return cls.get_filekey_from_pattern(key_string)

    @classmethod
    def get_filekey_from_filename(cls, filename):
        return cls.get_filekey_from_pattern(filename, processing_pattern())

    @classmethod
    def get_filekey_from_pattern(cls, filename, pattern=None):
        if isinstance(pattern, Path):
            pattern = pattern.as_posix()
        filename = str(filename)
        key_pattern_rx = re.compile(
            regex_from_filepattern(cls.key_pattern if pattern is None else pattern)
        )

        if key_pattern_rx.match(filename) is None:
            return None
        else:
            d = key_pattern_rx.match(filename).groupdict()
            for entry in list(d):
                if entry not in cls._fields:
                    d.pop(entry)
            for wildcard in cls._fields:
                if wildcard not in d:
                    d[wildcard] = "*"
            return cls(**d)

    @classmethod
    def unix_time_from_string(cls, value):
        key_class = cls.from_string(value)
        return unix_time(key_class.timestamp)

    def get_unix_timestamp(self):
        return unix_time(self.timestamp)

    @classmethod
    def parse_keypart(cls, keypart):
        keypart_rx = re.compile(cls.re_pattern)
        d = keypart_rx.match(keypart).groupdict()
        for key in d:
            if d[key] is None:
                d[key] = "*"
        return cls(**d)

    def expand(self, file_pattern, **kwargs):
        file_pattern = str(file_pattern)
        wildcard_dict = self._asdict()
        if kwargs is not None:
            for key, value in kwargs.items():
                wildcard_dict[key] = value
        wildcard_dict = {
            wildcard: [wildcard_value]
            if isinstance(wildcard_value, str)
            else wildcard_value
            for wildcard, wildcard_value in wildcard_dict.items()
        }
        formatter = string.Formatter()
        result = []
        for combo in product(*wildcard_dict.values()):
            substitution = dict(zip(list(wildcard_dict), combo))
            result.append(formatter.vformat(file_pattern, (), substitution))
        return result

    def get_path_from_filekey(self, pattern, **kwargs):
        if kwargs is None:
            return self.expand(pattern, **kwargs)
        else:
            for entry, value in kwargs.items():
                if isinstance(value, dict):
                    if len(next(iter(set(value).intersection(self._list())))) > 0:
                        kwargs[entry] = value[
                            next(iter(set(value).intersection(self._list())))
                        ]
                    else:
                        kwargs.pop(entry)
            return self.expand(pattern, **kwargs)

    # get_path_from_key
    @classmethod
    def get_full_path_from_filename(cls, filename, pattern, path_pattern):
        return cls.get_path_from_filekey(
            cls.get_filekey_from_pattern(filename, pattern), path_pattern
        )

    @staticmethod
    def tier_files(setup, keys, tier):
        fn_pattern = get_pattern_tier(setup, tier)
        files = []
        for line in keys:
            tier_filename = FileKey.get_full_path_from_filename(
                line, FileKey.key_pattern, fn_pattern
            )
            files += tier_filename
        return files


class ProcessingFileKey(FileKey):
    _fields = (*FileKey._fields, "processing_step")
    key_pattern = processing_pattern()

    def __new__(cls, experiment, period, run, datatype, timestamp, processing_step):
        self = super().__new__(cls, experiment, period, run, datatype, timestamp)
        if "_" in processing_step:
            splits = processing_step.split("_", 2)
            self.processing_type = splits[0]
            self.tier = splits[1]
            if len(splits) > 2:
                self.identifier = splits[2]
            else:
                self.identifier = None
        else:
            self.processing_type = processing_step
            self.tier = None
            self.identifier = None

        return self

    def _list(self):
        return [*super()._list(), self.processing_step]

    def _asdict(self):
        dic = super()._asdict()
        dic["processing_step"] = self.processing_step
        return dic

    @property
    def processing_step(self):
        if self.identifier is not None:
            return f"{self.processing_type}_{self.tier}_{self.identifier}"
        else:
            return f"{self.processing_type}_{self.tier}"

    @property
    def name(self):
        return f"{super().name}-{self.processing_step}"

    def get_path_from_filekey(self, pattern, **kwargs):
        if isinstance(pattern, Path):
            pattern = pattern.as_posix()
        if not isinstance(pattern, str):
            pattern = pattern(self.tier, self.identifier)
        if kwargs is None:
            return self.expand(pattern, **kwargs)
        else:
            for entry, value in kwargs.items():
                if isinstance(value, dict):
                    if len(next(iter(set(value).intersection(self._list())))) > 0:
                        kwargs[entry] = value[
                            next(iter(set(value).intersection(self._list())))
                        ]
                    else:
                        kwargs.pop(entry)
            return self.expand(pattern, **kwargs)


class ChannelProcKey(FileKey):
    re_pattern = "all(-(?P<experiment>[^-]+)(\\-(?P<period>[^-]+)(\\-(?P<run>[^-]+)(\\-(?P<datatype>[^-]+)(\\-(?P<timestamp>[^-]+)(\\-(?P<channel>[^-]+))?)?)?)?)?)?$"
    key_pattern = full_channel_pattern_with_extension()
    _fields = (*FileKey._fields, "channel")

    def __new__(cls, experiment, period, run, datatype, timestamp, channel):
        self = super().__new__(cls, experiment, period, run, datatype, timestamp)
        self.channel = channel
        return self

    @property
    def name(self):
        return f"{super().name}-{self.channel}"

    def _asdict(self):
        dic = super()._asdict()
        dic["channel"] = self.channel
        return dic

    @staticmethod
    def get_channel_files(keypart, par_pattern, chan_list):
        if isinstance(par_pattern, Path):
            par_pattern = par_pattern.as_posix()
        d = ChannelProcKey.parse_keypart(keypart)
        filenames = []
        for chan in chan_list:
            wildcards_dict = d._asdict()
            wildcards_dict.pop("channel")
            formatter = string.Formatter()
            wildcards_dict["channel"] = chan
            file = formatter.vformat(par_pattern, (), wildcards_dict)
            filenames.append(file)
        return filenames


def per_grouper(files):
    """
    Returns list containing lists of each run
    """

    pers = []
    per_files = []
    for file in files:
        fk = ProcessingFileKey.get_filekey_from_pattern(Path(file).name)
        if f"{fk.experiment}-{fk.period}" not in pers:
            pers.append(f"{fk.experiment}-{fk.period}")
            per_files.append([])
        for i, run in enumerate(pers):
            if run == f"{fk.experiment}-{fk.period}":
                per_files[i].append(file)
    return per_files


def run_grouper(files):
    """Returns list containing lists of each run"""
    runs = []
    run_files = []
    for file in files:
        fk = ProcessingFileKey.get_filekey_from_pattern(Path(file).name)
        if f"{fk.experiment}-{fk.period}-{fk.run}" not in runs:
            runs.append(f"{fk.experiment}-{fk.period}-{fk.run}")
            run_files.append([])
        for i, run in enumerate(runs):
            if run == f"{fk.experiment}-{fk.period}-{fk.run}":
                run_files[i].append(file)
    return run_files


def run_splitter(files):
    """
    Returns list containing lists of each run
    """

    runs = []
    run_files = []
    for file in files:
        fk = ProcessingFileKey.get_filekey_from_pattern(Path(file).name)
        if f"{fk.period}-{fk.run}" not in runs:
            runs.append(f"{fk.period}-{fk.run}")
            run_files.append([])
        for i, run in enumerate(runs):
            if run == f"{fk.period}-{fk.run}":
                run_files[i].append(file)
    return run_files
