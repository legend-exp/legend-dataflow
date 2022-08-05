import sys
import snakemake as smk
import re
from .patterns import *
from .utils import *
from collections import namedtuple


class FileKey(namedtuple('FileKey', ['experiment', 'period', 'run', 'datatype', 'timestamp'])):
    __slots__ = ()
    
    key_pattern = '(-(?P<experiment>[^-]+)(\\-(?P<period>[^-]+)(\\-(?P<run>[^-]+)(\\-(?P<datatype>[^-]+)(\\-(?P<timestamp>[^-]+))?)?)?)?)?$'
    
    @property
    def name(self):
        return f'{self.experiment}-{self.period}-{self.run}-{self.datatype}-{self.timestamp}'
        
    
    @property
    def key(self):
        return f'-{self.experiment}-{self.period}-{self.run}-{self.datatype}-{self.timestamp}'
    
    @classmethod
    def from_string(cls, key_string):
        parts = key_string.split("-")
        experiment, period, run, datatype, timestamp = parts
        return cls(experiment, period, run, datatype, timestamp)
    
    @classmethod
    def get_filekey_from_filename(cls, filename):
        experiment, period, run, datatype, timestamp, processing_step = filename.split("-")
        return cls(experiment, period, run, datatype, timestamp)
    
    @classmethod
    def get_filekey_from_pattern(cls, filename, pattern=None):
        if pattern is None:
            key_pattern_rx = re.compile(smk.io.regex(key_pattern()))
        else:
            key_pattern_rx = re.compile(smk.io.regex(pattern))
        if key_pattern_rx.match(filename) is None:
            print("pattern does not match")
        else:
            d = key_pattern_rx.match(filename).groupdict()
            for entry in list(d):
                if entry not in ['experiment', 'period', 'run', 'datatype', 'timestamp']:
                    d.pop(entry)
            for wildcard in ['experiment', 'period', 'run', 'datatype', 'timestamp']:
                if wildcard not in d:
                    d[wildcard] = '*'
            return cls(**d)
        
    def get_path_from_filekey(self, path_patttern):
        return smk.io.expand(path_patttern, **self._asdict())
    
    @staticmethod
    def get_full_path_from_filename(filename, pattern, path_pattern):
        return FileKey.get_path_from_filekey(FileKey.get_filekey_from_pattern(filename, pattern) , path_pattern)
    
    @classmethod
    def unix_time_from_string(cls, value):
        key_class = cls.from_string(value)
        return unix_time(key_class.timestamp)
        
    
    def get_unix_timestamp(self):
        return unix_time(self.timestamp)
    
    
    def get_file(self, pattern,name=None):
        if name is None:
            return smk.io.expand(pattern, **self._asdict())
        else:
            return smk.io.expand(pattern, **self._asdict(), name=name[self.datatype])
    
    @classmethod
    def parse_keypart(cls,keypart):
        keypart_rx = re.compile(FileKey.key_pattern)
        d = keypart_rx.match(keypart).groupdict()
        for key in d:
            if d[key] is None:
                d[key] = "*"
        return cls(**d)
    
    def __str__(self):
        return self.name

    @staticmethod
    def tier_files(setup, dataset_file, tier):
        fn_pattern = get_pattern_tier(setup, tier)
        files = []
        with open(dataset_file) as f:
            for line in f:
                tier_filename = FileKey.get_full_path_from_filename(line, key_pattern(), fn_pattern)
                files += tier_filename
        return files

class ProcessingFileKey(FileKey):
    #__slots__ = ()
    
    def __new__(cls, experiment, period, run, datatype, timestamp, processing_step):
        self = super(ProcessingFileKey,cls).__new__(cls, experiment, period, run, datatype, timestamp)
        self.processing_step = processing_step
        return self
    
    general_pattern = '{experiment}-{period}-{run}-{datatype}-{timestamp}-{processing_step}.{ext}'
    
    @property
    def name(self):
        return f'{self.experiment}-{self.period}-{self.run}-{self.datatype}-{self.timestamp}-{self.processing_step}'

    
    @classmethod
    def from_string(cls, key_string):
        parts = key_string.split("-")
        experiment, period, run, datatype, timestamp, processing_step = parts
        return cls(experiment, period, run, datatype, timestamp, processing_step)
    
    @classmethod
    def get_filekey_from_filename(cls, filename):
        experiment, period, run, datatype, timestamp, processing_step = filename.split("-")
        return cls(experiment, period, run, datatype, timestamp, processing_step.split(".")[0])
    
    @classmethod
    def get_filekey_from_pattern(cls, filename, pattern=None):
        if pattern is None:
            key_pattern_rx = re.compile(smk.io.regex(ProcessingFileKey.general_pattern))
        else:
            key_pattern_rx = re.compile(smk.io.regex(pattern))
        if key_pattern_rx.match(filename) is None:
            print("pattern does not match")
        else:
            d = key_pattern_rx.match(filename).groupdict()
            for entry in list(d):
                if entry not in ['experiment', 'period', 'run', 'datatype', 'timestamp', 'processing_step']:
                    d.pop(entry)
            return cls(**d)    
        
class ChannelFileKey(FileKey):
    
    channel_pattern = 'all(-(?P<experiment>[^-]+)(\\-(?P<period>[^-]+)(\\-(?P<run>[^-]+)(\\-(?P<datatype>[^-]+)(\\-(?P<timestamp>[^-]+)(\\-(?P<channel>[^-]+))?)?)?)?)?)?$'
    
    def __new__(cls, experiment, period, run, datatype,timestamp, channel):
        self = super(ChannelFileKey,cls).__new__(cls, experiment, period, run, datatype,timestamp)
        self.channel = channel
        return self
    
    @property
    def name(self):
        return f'{self.experiment}-{self.period}-{self.run}-{self.datatype}-{self.timestamp}-{self.channel}'
    
    @classmethod
    def parse_keypart(cls,keypart):
        keypart_rx = re.compile(ChannelFileKey.channel_pattern)
        d = keypart_rx.match(keypart).groupdict()
        for key in d:
            if d[key] is None:
                d[key] = "*"
        return cls(**d)
    
    def __dict__(self):
        out = self._asdict()
        out["channel"] = self.channel
        return out
    
    @staticmethod
    def get_channel_files(setup, keypart, tier, dataset_file):
        d = ChannelFileKey.parse_keypart(keypart)
        par_pattern = get_pattern_pars_tmp_channel(setup,tier)
        filenames = []
        with open(dataset_file) as f:
            chan_list = f.read().splitlines()
            for chan in chan_list:
                wildcards_dict = d.__dict__()
                wildcards_dict.pop("channel")
                file = smk.io.expand(par_pattern, **wildcards_dict , channel = chan)[0]
                filenames.append(file)
        return filenames
    
    @classmethod
    def get_filekey_from_pattern(cls, filename, pattern=None):
        if pattern is None:
            key_pattern_rx = re.compile(smk.io.regex(key_pattern()))
        else:
            key_pattern_rx = re.compile(smk.io.regex(pattern))
        if key_pattern_rx.match(filename) is None:
            print("pattern does not match")
        else:
            d = key_pattern_rx.match(filename).groupdict()
            for entry in list(d):
                if entry not in ['experiment', 'period', 'run', 'datatype', 'timestamp', 'channel']:
                    d.pop(entry)
            for wildcard in ['experiment', 'period', 'run', 'datatype', 'timestamp', 'channel']:
                if wildcard not in d:
                    d[wildcard] = '*'
            return cls(**d)

class ChannelProcessingFileKey(FileKey):
    
    def __new__(cls, experiment, period, run, datatype,timestamp, channel, processing_step):
        self = super(ChannelProcessingFileKey,cls).__new__(cls, experiment, period, run, datatype,timestamp)
        self.channel = channel
        self.processing_step = processing_step
        return self
    
    @property
    def name(self):
        return f'{self.experiment}-{self.period}-{self.run}-{self.datatype}-{self.timestamp}-{self.channel}-{self.processing_step}'
    
    @classmethod
    def get_filekey_from_pattern(cls, filename, pattern=None):
        if pattern is None:
            key_pattern_rx = re.compile(smk.io.regex(full_channel_pattern()))
        else:
            key_pattern_rx = re.compile(smk.io.regex(pattern))
        if key_pattern_rx.match(filename) is None:
            print("pattern does not match")
        else:
            d = key_pattern_rx.match(filename).groupdict()
            for entry in list(d):
                if entry not in ['experiment', 'period', 'run', 'datatype', 'timestamp', 'channel','processing_step']:
                    d.pop(entry)
            for wildcard in ['experiment', 'period', 'run', 'datatype', 'timestamp', 'channel','processing_step']:
                if wildcard not in d:
                    d[wildcard] = '*'
            return cls(**d)   
    
    def get_path_from_filekey(self, path_patttern):
        return smk.io.expand(path_patttern, **self._asdict(), **self.__dict__)