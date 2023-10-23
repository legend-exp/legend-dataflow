import re
from collections import namedtuple

import snakemake as smk

from .patterns import *
from .utils import *

# key_pattern -> key
#


class FileKey(namedtuple("FileKey", ["experiment", "period", "run", "datatype", "timestamp"])):
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
        if pattern is None:
            key_pattern_rx = re.compile(smk.io.regex(cls.key_pattern))
        else:
            key_pattern_rx = re.compile(smk.io.regex(pattern))
        if key_pattern_rx.match(filename) is None:
            pass
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

    def get_path_from_filekey(self, pattern, **kwargs):
        if kwargs is None:
            return smk.io.expand(pattern, **self._asdict())
        else:
            for entry, value in kwargs.items():
                if isinstance(value, dict):
                    if len(list(set(value).intersection(self._list()))[0]) > 0:
                        kwargs[entry] = value[list(set(value).intersection(self._list()))[0]]
                    else:
                        kwargs.pop(entry)
            return smk.io.expand(pattern, **self._asdict(), **kwargs)

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
        if not isinstance(pattern, str):
            pattern = pattern(self.tier, self.identifier)
        if kwargs is None:
            return smk.io.expand(pattern, **self._asdict())
        else:
            for entry, value in kwargs.items():
                if isinstance(value, dict):
                    if len(list(set(value).intersection(self._list()))[0]) > 0:
                        kwargs[entry] = value[list(set(value).intersection(self._list()))[0]]
                    else:
                        kwargs.pop(entry)
            return smk.io.expand(pattern, **self._asdict(), **kwargs)


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
    def get_channel_files(setup, keypart, tier, chan_list, name=None):
        d = ChannelProcKey.parse_keypart(keypart)
        par_pattern = get_pattern_pars_tmp_channel(setup, tier, name)
        filenames = []
        for chan in chan_list:
            wildcards_dict = d._asdict()
            wildcards_dict.pop("channel")
            file = smk.io.expand(par_pattern, **wildcards_dict, channel=chan)[0]
            filenames.append(file)
        return filenames
