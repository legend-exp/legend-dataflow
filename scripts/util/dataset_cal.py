"""
This module uses the partition database files to the necessary inputs for partition calibrations
"""

import json
import os

from .CalibCatalog import PropsStream
from .FileKey import ChannelProcKey, ProcessingFileKey
from .patterns import (
    get_pattern_log_channel,
    get_pattern_pars_tmp_channel,
    get_pattern_plts_tmp_channel,
)
from .utils import filelist_path


class dataset_file:
    def __init__(self, setup, input_file):
        with open(input_file) as r:
            self.datasets = json.load(r)
        self.setup = setup

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
                    os.path.join(
                        filelist_path(self.setup),
                        f"all-{experiment}-{per}-*-{datatype}-{tier}.filelist",
                    )
                ]
            else:
                files += [
                    os.path.join(
                        filelist_path(self.setup),
                        f"all-{experiment}-{per}-{run}-{datatype}-{tier}.filelist",
                    )
                    for run in dataset[per]
                ]
        return files

    def get_par_files(
        self,
        catalog_file,
        dataset,
        channel,
        tier,
        experiment="l200",
        datatype="cal",
        name=None,
        extension="json",
    ):
        dataset = self.get_dataset(dataset, channel)
        all_par_files = []
        for item in PropsStream.read_from(catalog_file):
            par_files = item["apply"]
            for par_file in par_files:
                if par_file.split("-")[-1] == f"par_{tier}.json":
                    all_par_files.append(par_file)
        if channel == "default":
            channel = "{channel}"
        selected_par_files = []
        for par_file in all_par_files:
            fk = ProcessingFileKey.get_filekey_from_pattern(os.path.basename(par_file))
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
        catalog_file,
        dataset,
        channel,
        tier,
        experiment="l200",
        datatype="cal",
        name=None,
    ):
        dataset = self.get_dataset(dataset, channel)
        all_par_files = []
        for item in PropsStream.read_from(catalog_file):
            par_files = item["apply"]
            for par_file in par_files:
                if par_file.split("-")[-1] == f"par_{tier}.json":
                    all_par_files.append(par_file)
        if channel == "default":
            channel = "{channel}"
        selected_par_files = []
        for par_file in all_par_files:
            fk = ProcessingFileKey.get_filekey_from_pattern(os.path.basename(par_file))
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
        catalog_file,
        dataset,
        channel,
        tier,
        experiment="l200",
        datatype="cal",
        name=None,
    ):
        par_files = self.get_par_files(
            catalog_file,
            dataset,
            channel,
            tier,
            experiment=experiment,
            datatype=datatype,
            name=name,
        )
        fk = ChannelProcKey.get_filekey_from_pattern(os.path.basename(par_files[0]))
        if channel == "default":
            fk.channel = "{channel}"
        else:
            fk.channel = channel
        return fk.get_path_from_filekey(get_pattern_log_channel(self.setup, name))[0]

    def get_timestamp(
        self, catalog_file, dataset, channel, tier, experiment="l200", datatype="cal"
    ):
        par_files = self.get_par_files(
            catalog_file,
            dataset,
            channel,
            tier,
            experiment=experiment,
            datatype=datatype,
            name=None,
        )
        fk = ChannelProcKey.get_filekey_from_pattern(os.path.basename(par_files[0]))
        return fk.timestamp
