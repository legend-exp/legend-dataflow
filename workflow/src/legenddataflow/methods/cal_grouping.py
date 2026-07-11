"""
This module uses the partition database files to the necessary inputs for partition calibrations
"""

from __future__ import annotations

import logging
import re
from pathlib import Path

from dbetto import Props

from .FileKey import ChannelProcKey, ProcessingFileKey
from .paths import filelist_path
from .patterns import (
    get_pattern_log_channel,
    get_pattern_pars,
    get_pattern_pars_tmp_channel,
    get_pattern_plts_tmp_channel,
)

log = logging.getLogger(__name__)


class CalGrouping:
    """Resolve partition-level calibration groupings (from
    ``cal_groupings.yaml``): which runs belong to each partition per channel,
    and the corresponding filelists, parameter/plot files, log paths and
    wildcard constraints."""

    def __init__(self, setup, input_file):
        self.datasets = Props.read_from(input_file)
        self.expand_runs()
        self.setup = setup

    @staticmethod
    def _expand_range(runs, channel, part, per):
        if re.fullmatch(r"r\d+\.\.r\d+", runs) is None:
            msg = (
                f"malformed run range {runs!r} for channel {channel!r}, "
                f"partition {part!r}, period {per!r}: expected 'rNNN..rNNN'"
            )
            raise ValueError(msg)
        start, end = runs.split("..")
        return [f"r{x:03}" for x in range(int(start[1:]), int(end[1:]) + 1)]

    def expand_runs(self):
        """Expand ``r000..r005`` run-range strings into explicit run lists."""
        for channel, chan_dict in self.datasets.items():
            for part, part_dict in chan_dict.items():
                for per, runs in part_dict.items():
                    if isinstance(runs, str) and ".." in runs:
                        self.datasets[channel][part][per] = self._expand_range(
                            runs, channel, part, per
                        )
                    if isinstance(runs, list):
                        final_runs = []
                        for run in runs:
                            if ".." in run:
                                final_runs += self._expand_range(
                                    run, channel, part, per
                                )
                            else:
                                if re.fullmatch(r"r\d+", run) is None:
                                    log.warning(
                                        "run %r for channel %r, partition %r, "
                                        "period %r does not look like 'rNNN'",
                                        run,
                                        channel,
                                        part,
                                        per,
                                    )
                                final_runs.append(run)
                        self.datasets[channel][part][per] = final_runs

    def get_dataset(self, dataset, channel):
        """Return the ``{period: runs}`` dict for a partition, with
        channel-specific overrides applied on top of ``default``."""
        if "default" not in self.datasets:
            msg = "cal-groupings config has no 'default' section"
            raise KeyError(msg)
        partition_dict = self.datasets["default"].copy()
        if channel in self.datasets:
            partition_dict.update(self.datasets[channel])
        if dataset not in partition_dict:
            msg = (
                f"partition {dataset!r} not defined for channel {channel!r} "
                f"(or 'default'): available partitions are {sorted(partition_dict)}"
            )
            raise KeyError(msg)
        return partition_dict[dataset]

    def get_filelists(self, dataset, channel, tier, experiment="l200", datatype="cal"):
        """Return the filelist paths covering all runs in the partition."""
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
        pattern_func=get_pattern_pars_tmp_channel,
    ):
        """Select from ``catalog`` the per-channel parameter files whose keys
        fall inside the partition's periods and runs, expanded with
        ``pattern_func``."""
        dataset = self.get_dataset(dataset, channel)
        all_par_files = []
        for item in catalog.entries["all"]:
            par_files = item.file
            for par_file in par_files:
                if (
                    par_file.split("-")[-1]
                    == str(
                        get_pattern_pars(self.setup, tier, check_in_cycle=False).name
                    ).split("-")[-1]
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
                            pattern_func(
                                self.setup, tier, name=name, extension=extension
                            ),
                            channel=channel,
                        )[0]
                    )
                else:
                    selected_par_files.append(
                        fk.get_path_from_filekey(
                            pattern_func(
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
        """Like :meth:`get_par_files` but for the pickled plot outputs."""
        return self.get_par_files(
            catalog,
            dataset,
            channel,
            tier,
            experiment=experiment,
            datatype=datatype,
            name=name,
            extension="pkl",
            pattern_func=get_pattern_plts_tmp_channel,
        )

    def get_log_file(
        self,
        catalog,
        dataset,
        channel,
        tier,
        processing_timestamp,
        experiment="l200",
        datatype="cal",
        name=None,
    ):
        """Build the log-file path for the partition job, derived from the
        first matching parameter file."""
        par_files = self.get_par_files(
            catalog,
            dataset,
            channel,
            tier,
            experiment=experiment,
            datatype=datatype,
            name=name,
        )
        if len(par_files) > 0:
            fk = ChannelProcKey.get_filekey_from_pattern(Path(par_files[0]).name)
            if channel == "default":
                fk.channel = "{channel}"
            else:
                fk.channel = channel
            return fk.get_path_from_filekey(
                get_pattern_log_channel(self.setup, name, processing_timestamp)
            )[0]
        log.warning(
            "partition %r channel %r tier %r resolved to no par files; "
            "using placeholder log file name",
            dataset,
            channel,
            tier,
        )
        return "log.log"

    def get_timestamp(
        self, catalog, dataset, channel, tier, experiment="l200", datatype="cal"
    ):
        """Return the timestamp of the partition's first parameter file."""
        par_files = self.get_par_files(
            catalog,
            dataset,
            channel,
            tier,
            experiment=experiment,
            datatype=datatype,
            name=None,
        )
        if len(par_files) > 0:
            fk = ChannelProcKey.get_filekey_from_pattern(Path(par_files[0]).name)
            return fk.timestamp
        log.warning(
            "partition %r channel %r tier %r resolved to no par files; "
            "using placeholder timestamp",
            dataset,
            channel,
            tier,
        )
        return "20000101T000000Z"

    def get_wildcard_constraints(self, dataset, channel):
        """Build the channel wildcard regex; for ``default`` it excludes
        channels that have their own override for the partition's runs."""
        if channel == "default":
            exclude_chans = []
            default_runs = self.get_dataset(dataset, channel)
            for chan, chan_dict in self.datasets.items():
                if chan != "default":
                    for _, dataset_dict in chan_dict.items():
                        for period, runs in dataset_dict.items():
                            if period in default_runs:
                                for run in runs:
                                    if run in default_runs[period]:
                                        exclude_chans.append(chan)
            exclude_chans = list(set(exclude_chans))
            out_string = r"(?:"
            if len(exclude_chans) > 0:
                out_string += r"(?!(?:"
                out_string += exclude_chans[0]
                for chan in exclude_chans[1:]:
                    out_string += f"|{chan}"
                out_string += r"))"
            return out_string + r"[PCVB]\d\w{5})"
        return r"[PCVB]{1}\d{1}\w{5}"
