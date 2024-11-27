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

"""
This module stores the scripts for leading validity files based on timestamp and system
"""

import bisect
import collections
import copy
import json
import types
from collections import namedtuple
from pathlib import Path

import yaml

from .utils import unix_time


class Props:
    @staticmethod
    def read_from(sources):
        def read_impl(sources):
            if isinstance(sources, (str, Path)):
                file_name = sources
                if isinstance(file_name, str):
                    file_name = Path(file_name)
                if file_name.suffix in (".yaml", ".yml"):
                    with file_name.open() as file:
                        return yaml.safe_load(file)
                elif file_name.suffix == ".json":
                    with open(file_name) as file:
                        return json.load(file)
                else:
                    msg = f"Can't run Props.read_from on file with suffix {file_name.suffix}"
                    raise ValueError(msg)
            elif isinstance(sources, list):
                result = {}
                for p in map(read_impl, sources):
                    Props.add_to(result, p)
                return result
            else:
                msg = f"Can't run Props.read_from on sources-value of type {type(sources)}"
                raise ValueError(msg)

        return read_impl(sources)

    @staticmethod
    def add_to(props_a, props_b):
        a = props_a
        b = props_b

        for key in b:
            if key in a:
                if isinstance(a[key], dict) and isinstance(b[key], dict):
                    Props.add_to(a[key], b[key])
                elif a[key] != b[key]:
                    a[key] = copy.copy(b[key])
            else:
                a[key] = copy.copy(b[key])


class PropsStream:
    """Simple class to control loading of validity.yaml files"""

    @staticmethod
    def get(value):
        if isinstance(value, str):
            return PropsStream.read_from(value)

        if isinstance(value, (collections.abc.Sequence, types.GeneratorType)):
            return value

        msg = f"Can't get PropsStream from value of type {type(value)}"
        raise ValueError(msg)

    @staticmethod
    def read_from(file_name):
        with Path(file_name).open() as r:
            file = yaml.safe_load(r)
        file = sorted(file, key=lambda item: unix_time(item["valid_from"]))
        yield from file


class Catalog(namedtuple("Catalog", ["entries"])):
    """Implementation of the `YAML metadata validity specification <https://legend-exp.github.io/legend-data-format-specs/dev/metadata/#Specifying-metadata-validity-in-time-(and-system)>`_."""

    __slots__ = ()

    class Entry(namedtuple("Entry", ["valid_from", "file"])):
        __slots__ = ()

    @staticmethod
    def get(value):
        if isinstance(value, Catalog):
            return value

        if isinstance(value, str):
            return Catalog.read_from(value)

        msg = f"Can't get Catalog from value of type {type(value)}"
        raise ValueError(msg)

    @staticmethod
    def read_from(file_name):
        """Read from a valdiity YAML file and build a Catalog object"""
        entries = {}
        for props in PropsStream.get(file_name):
            timestamp = props["valid_from"]
            system = "all" if props.get("category") is None else props["category"]
            file_key = props["apply"]
            if system not in entries:
                entries[system] = []
            mode = "append" if props.get("mode") is None else props["mode"]
            mode = "reset" if len(entries[system]) == 0 else mode
            if mode == "reset":
                new = file_key
            elif mode == "append":
                new = entries[system][-1].file.copy() + file_key
            elif mode == "remove":
                new = entries[system][-1].file.copy()
                for file in file_key:
                    new.remove(file)
            elif mode == "replace":
                new = entries[system][-1].file.copy()
                if len(file_key) != 2:
                    msg = f"Invalid number of elements in replace mode: {len(file_key)}"
                    raise ValueError(msg)
                new.remove(file_key[0])
                new += [file_key[1]]

            else:
                msg = f"Unknown mode for {timestamp}"
                raise ValueError(msg)

            if timestamp in [entry.valid_from for entry in entries[system]]:
                msg = (
                    f"Duplicate timestamp: {timestamp}, use reset mode instead with a single entry"
                )
                raise ValueError(msg)
            entries[system].append(Catalog.Entry(unix_time(timestamp), new))

        for system in entries:
            entries[system] = sorted(entries[system], key=lambda entry: entry.valid_from)
        return Catalog(entries)

    def valid_for(self, timestamp, system="all", allow_none=False):
        """Get the valid entries for a given timestamp and system"""
        if system in self.entries:
            valid_from = [entry.valid_from for entry in self.entries[system]]
            pos = bisect.bisect_right(valid_from, unix_time(timestamp))
            if pos > 0:
                return self.entries[system][pos - 1].file

            if system != "all":
                return self.valid_for(timestamp, system="all", allow_none=allow_none)

            if allow_none:
                return None

            msg = f"No valid entries found for timestamp: {timestamp}, system: {system}"
            raise RuntimeError(msg)

        if system != "all":
            return self.valid_for(timestamp, system="all", allow_none=allow_none)

        if allow_none:
            return None

        msg = f"No entries found for system: {system}"
        raise RuntimeError(msg)

    @staticmethod
    def get_files(catalog_file, timestamp, category="all"):
        """Helper function to get the files for a given timestamp and category"""
        catalog = Catalog.read_from(catalog_file)
        return Catalog.valid_for(catalog, timestamp, category)
