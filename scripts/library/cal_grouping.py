"""
This module uses the partition database files to the necessary inputs for partition calibrations
"""

import json
from pathlib import Path

import yaml

from .FileKey import ChannelProcKey, ProcessingFileKey
from .patterns import (
    get_pattern_log_channel,
    get_pattern_pars,
    get_pattern_pars_tmp_channel,
    get_pattern_plts_tmp_channel,
)
from .utils import filelist_path


class CalGrouping:
    def __init__(self, setup, input_file):
        if Path(input_file).suffix == ".json":
            with Path(input_file).open() as r:
                self.datasets = json.load(r)
        elif Path(input_file).suffix in (".yaml", ".yml"):
            with Path(input_file).open() as r:
                self.datasets = yaml.safe_load(r)
        self.expand_runs()
        self.setup = setup

    def expand_runs(self):
        for channel, chan_dict in self.datasets.items():
            for part, part_dict in chan_dict.items():
                for per, runs in part_dict.items():
                    if isinstance(runs, str) and ".." in runs:
                        start, end = runs.split("..")
                        self.datasets[channel][part][per] = [
                            f"r{x:03}" for x in range(int(start[1:]), int(end[1:]) + 1)
                        ]

    def get_dataset(self, dataset, channel):
        partition_dict = self.datasets["default"].copy()
        if channel in self.datasets:
            partition_dict.update(self.datasets[channel])
        return partition_dict[dataset]

    def get_filelists(self, dataset, channel, tier, experiment="l200", datatype="cal"):
        dataset = self.get_dataset(dataset, channel)
        files = []
        for per in dataset:
            if dataset[per] == "all":
                files += [
                    Path(filelist_path(self.setup))
                    / f"all-{experiment}-{per}-*-{datatype}-{tier}.filelist"
                ]
            else:
                files += [
                    Path(filelist_path(self.setup))
                    / f"all-{experiment}-{per}-{run}-{datatype}-{tier}.filelist"
                    for run in dataset[per]
                ]
        return files

    def get_par_files(
        self,
        catalog,
        dataset,
        channel,
        tier,
        experiment="l200",
        datatype="cal",
        name=None,
        extension="yaml",
    ):
        dataset = self.get_dataset(dataset, channel)
        all_par_files = []
        for item in catalog:
            par_files = item.apply
            for par_file in par_files:
                if (
                    par_file.split("-")[-1]
                    == str(get_pattern_pars(self.setup, tier, check_in_cycle=False).name).split(
                        "-"
                    )[-1]
                ):
                    all_par_files.append(par_file)
        if channel == "default":
            channel = "{channel}"
        selected_par_files = []
        for par_file in all_par_files:
            fk = ProcessingFileKey.get_filekey_from_pattern(Path(par_file).name)
            if (
                fk.datatype == datatype
                and fk.experiment == experiment
                and fk.period in list(dataset)
                and (dataset[fk.period] == "all" or fk.run in dataset[fk.period])
            ):
                if name is not None:
                    selected_par_files.append(
                        fk.get_path_from_filekey(
                            get_pattern_pars_tmp_channel(
                                self.setup, tier, name=name, extension=extension
                            ),
                            channel=channel,
                        )[0]
                    )
                else:
                    selected_par_files.append(
                        fk.get_path_from_filekey(
                            get_pattern_pars_tmp_channel(
                                self.setup, tier, name=name, extension=extension
                            ),
                            channel=channel,
                        )[0]
                    )
        return selected_par_files

    def get_plt_files(
        self,
        catalog,
        dataset,
        channel,
        tier,
        experiment="l200",
        datatype="cal",
        name=None,
    ):
        dataset = self.get_dataset(dataset, channel)
        all_par_files = []
        for item in catalog:
            par_files = item.apply
            for par_file in par_files:
                if (
                    par_file.split("-")[-1]
                    == str(get_pattern_pars(self.setup, tier, check_in_cycle=False).name).split(
                        "-"
                    )[-1]
                ):
                    all_par_files.append(par_file)
        if channel == "default":
            channel = "{channel}"
        selected_par_files = []
        for par_file in all_par_files:
            fk = ProcessingFileKey.get_filekey_from_pattern(Path(par_file).name)
            if (
                fk.datatype == datatype
                and fk.experiment == experiment
                and fk.period in list(dataset)
                and (dataset[fk.period] == "all" or fk.run in dataset[fk.period])
            ):
                if name is not None:
                    selected_par_files.append(
                        fk.get_path_from_filekey(
                            get_pattern_plts_tmp_channel(self.setup, tier, name=name),
                            channel=channel,
                        )[0]
                    )
                else:
                    selected_par_files.append(
                        fk.get_path_from_filekey(
                            get_pattern_plts_tmp_channel(self.setup, tier, name=name),
                            channel=channel,
                        )[0]
                    )
        return selected_par_files

    def get_log_file(
        self,
        catalog,
        dataset,
        channel,
        tier,
        experiment="l200",
        datatype="cal",
        name=None,
    ):
        par_files = self.get_par_files(
            catalog,
            dataset,
            channel,
            tier,
            experiment=experiment,
            datatype=datatype,
            name=name,
        )
        fk = ChannelProcKey.get_filekey_from_pattern(Path(par_files[0]).name)
        if channel == "default":
            fk.channel = "{channel}"
        else:
            fk.channel = channel
        return fk.get_path_from_filekey(get_pattern_log_channel(self.setup, name))[0]

    def get_timestamp(self, catalog, dataset, channel, tier, experiment="l200", datatype="cal"):
        par_files = self.get_par_files(
            catalog,
            dataset,
            channel,
            tier,
            experiment=experiment,
            datatype=datatype,
            name=None,
        )
        fk = ChannelProcKey.get_filekey_from_pattern(Path(par_files[0]).name)
        return fk.timestamp

    def get_wildcard_constraints(self, dataset, channel):
        if channel == "default":
            exclude_chans = []
            default_runs = self.get_dataset(dataset, channel)
            for channel, chan_dict in self.datasets.items():
                if channel != "default":
                    for _, dataset_dict in chan_dict.items():
                        for period, runs in dataset_dict.items():
                            if period in default_runs:
                                for run in runs:
                                    if run in default_runs[period]:
                                        exclude_chans.append(channel)
            exclude_chans = set(exclude_chans)
            out_string = ""
            for channel in exclude_chans:
                out_string += f"(?!{channel})"
            return out_string + r"^[VPCB]\d{1}\w{5}$"
        else:
            return r"^[VPCB]\d{1}\w{5}$"
