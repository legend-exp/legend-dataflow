import sys
import snakemake as smk
import re
from .patterns import *
from .utils import *
from collections import namedtuple

# key_pattern -> key
# 


class FileKey(namedtuple('FileKey', ['experiment', 'period', 'run', 'datatype', 'timestamp'])):
    __slots__ = ()
    
    re_pattern = '(-(?P<experiment>[^-]+)(\\-(?P<period>[^-]+)(\\-(?P<run>[^-]+)(\\-(?P<datatype>[^-]+)(\\-(?P<timestamp>[^-]+))?)?)?)?)?$'
    key_pattern = key_pattern()
    
    @property
    def name(self):
        return f'{self.experiment}-{self.period}-{self.run}-{self.datatype}-{self.timestamp}'
        
    @property
    def key(self):
        return f'-{self.name}'
        
    def __str__(self):
        return self.name
    
    @classmethod
    def from_string(cls, key_string):
        return cls.get_filekey_from_pattern(key_string)
    
    @classmethod
    def get_filekey_from_filename(cls, filename):
        return cls.get_filekey_from_pattern(key_string, processing_pattern())
    
    @classmethod
    def get_filekey_from_pattern(cls, filename, pattern=None):
        if pattern is None:
            key_pattern_rx = re.compile(smk.io.regex(cls.key_pattern))
        else:
            key_pattern_rx = re.compile(smk.io.regex(pattern))
        if key_pattern_rx.match(filename) is None:
            print("pattern does not match")
        else:
            d = key_pattern_rx.match(filename).groupdict()
            for entry in list(d):
                if entry not in cls._fields:
                    d.pop(entry)
            for wildcard in cls._fields:
                if wildcard not in d:
                    d[wildcard] = '*'
            return cls(**d)
        
    @classmethod
    def unix_time_from_string(cls, value):
        key_class = cls.from_string(value)
        return unix_time(key_class.timestamp)
        
    
    def get_unix_timestamp(self):
        return unix_time(self.timestamp)    
    
    @classmethod
    def parse_keypart(cls,keypart):
        keypart_rx = re.compile(cls.re_pattern)
        d = keypart_rx.match(keypart).groupdict()
        for key in d:
            if d[key] is None:
                d[key] = "*"
        return cls(**d)
        
    def get_path_from_filekey(self, pattern, name=None):
        if name is None:
            return smk.io.expand(pattern, **self._asdict())
        else:
            return smk.io.expand(pattern, **self._asdict(), name=name[self.datatype])
    
    # get_path_from_key
    @classmethod
    def get_full_path_from_filename(cls, filename, pattern, path_pattern):
        return cls.get_path_from_filekey(cls.get_filekey_from_pattern(filename, pattern) , path_pattern)

    @staticmethod
    def tier_files(setup, dataset_file, tier):
        fn_pattern = get_pattern_tier(setup, tier)
        files = []
        with open(dataset_file) as f:
            for line in f:
                tier_filename = FileKey.get_full_path_from_filename(line, FileKey.key_pattern, fn_pattern)
                files += tier_filename
        return files
    
class ProcessingFileKey(FileKey):
    _fields = FileKey._fields+("processing_step", )
    pattern = processing_pattern()
    
    def __new__(cls, experiment, period, run, datatype, timestamp, processing_step):
        self = super(ProcessingFileKey, cls).__new__(cls, experiment, period, run, datatype, timestamp)
        if "_" in processing_step:
            splits = processing_step.split("_",2)
            self.processing_type = splits[0]
            self.tier = splits[1]
            if len(splits)>2:
                self.identifier = splits[2]
            else:
                self.identifier=None
        else:
            self.processing_type = processing_step
            self.tier=None
            self.identifier=None
        
        return self
    
    general_pattern = processing_pattern()
    
    @property
    def processing_step(self):
        if self.identifier is not None:
            return f"{self.processing_type}_{self.tier}_{self.identifier}"
        else:
            return f"{self.processing_type}_{self.tier}"
    
    @property
    def name(self):
        if self.identifier is not None:
            return f"{super().name}-{self.processing_step}"
        else:
            return f"{super().name}-{self.processing_step}"
    
    def get_path_from_filekey(self, pattern, name=None):
        if isinstance(pattern, str):
            if name is None:
                return smk.io.expand(pattern, **self._asdict())
            else:
                return smk.io.expand(pattern, **self._asdict(), name=name[self.datatype])
        else:
            final_pattern = pattern(self.tier, self.identifier)
            return smk.io.expand(final_pattern, **self._asdict())

    
class ChannelProcKey(ProcessingFileKey):
    
    _fields = ProcessingFileKey._fields+("channel",)
    
    def __new__(cls, experiment, period, run, datatype,timestamp, channel, processing_step):
        self = super(ChannelProcKey,cls).__new__(cls, experiment, period, run, datatype,timestamp, processing_step)
        self.channel = channel
        return self
    
    @property
    def name(self):
        return f'{super(Filekey).name}-{self.channel}-{super().processing_step}'