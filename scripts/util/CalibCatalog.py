# -*- coding: utf-8 -*-
#
# Copyright (C) 2015 Oliver Schulz <oschulz@mpp.mpg.de>
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
# http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#

from collections import namedtuple
import bisect
import types
import collections
import json
import copy
import os
from string import Template
from .utils import *

class PropsStream():
    @staticmethod
    def get(value):
        if isinstance(value, str):
            return PropsStream.read_from(value)
        elif isinstance(value, collections.Sequence) or isinstance(value, types.GeneratorType):
            return value
        else:
            raise ValueError("Can't get PropsStream from value of type {t}".format(t = type(source)))


    @staticmethod
    def read_from(file_name):
        with open(file_name, 'r') as file:
            for json_str in file:
                yield json.loads(json_str)

class CalibCatalog(namedtuple('CalibCatalog', ['entries'])):
    __slots__ = ()

    class Entry(namedtuple('Entry', ['valid_from','file'])):
        __slots__ = ()

    @staticmethod
    def read_from(file_name):
        entries = {}

        for props in PropsStream.get(file_name):
            timestamp = props["valid_from"]
            if props.get("select") is None:
                system = "all"
            else:
                system = props["select"]
            file_key = props["apply"]
            if system not in entries:
                entries[system] = []
            entries[system].append(CalibCatalog.Entry(unix_time(timestamp),file_key))

        for system in entries:
            entries[system] = sorted(
                entries[system],
                key = lambda entry: entry.valid_from
            )
        return CalibCatalog(entries)


    def calib_for(self, timestamp, system="all", allow_none = False):
        if system in self.entries:
            valid_from = [ entry.valid_from for entry in self.entries[system]]
            pos = bisect.bisect_right(valid_from, unix_time(timestamp))
            if pos > 0:
                return self.entries[system][pos - 1].file
            else:
                if allow_none: return None
                else: raise RuntimeError(f'No valid calibration found for timestamp: {timestamp}, system: {system}')
        else:
            if allow_none: return None
            else: raise RuntimeError(f'No calibrations found for system: {system}')
    
    @staticmethod
    def get_calib_files(catalog_file, timestamp, category="all"):
        catalog = CalibCatalog.read_from(catalog_file)
        return CalibCatalog.calib_for(catalog,timestamp, category)
    