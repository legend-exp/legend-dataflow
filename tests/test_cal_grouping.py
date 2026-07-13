from __future__ import annotations

import json
import logging
from pathlib import Path

import pytest
import yaml
from dbetto.catalog import Catalog
from legenddataflowscripts.workflow import subst_vars

from legenddataflow.methods import CalGrouping

testprod = Path(__file__).parent / "dummy_cycle"

with (testprod / "config.yaml").open() as r:
    setup = yaml.safe_load(r)
subst_vars(setup, var_values={"_": str(testprod)})


catalog = Catalog.get(
    [
        {
            "valid_from": "20210101T000000Z",
            "apply": ["l200-p01-r001-cal-20210101T000000Z-par_dsp.yaml"],
        },
        {
            "valid_from": "20210102T000000Z",
            "apply": ["l200-p01-r002-cal-20210102T000000Z-par_dsp.yaml"],
            "mode": "reset",
        },
    ]
)


@pytest.fixture
def input_file_json(tmp_path):
    data = {"default": {"calgroup001a": {"p01": "r001..r003"}}}
    file = tmp_path / "input.json"
    file.write_text(json.dumps(data))
    return file


@pytest.fixture
def input_file_yaml(tmp_path):
    data = {"default": {"calgroup001a": {"p01": "r001..r003"}}}
    file = tmp_path / "input.yaml"
    file.write_text(yaml.dump(data))
    return file


def test_expand_runs_json(input_file_json):
    cal_grouping = CalGrouping(setup, input_file_json)
    assert cal_grouping.datasets["default"]["calgroup001a"]["p01"] == [
        "r001",
        "r002",
        "r003",
    ]


def test_expand_runs_yaml(input_file_yaml):
    cal_grouping = CalGrouping(setup, input_file_yaml)
    assert cal_grouping.datasets["default"]["calgroup001a"]["p01"] == [
        "r001",
        "r002",
        "r003",
    ]


def test_get_dataset(input_file_json):
    cal_grouping = CalGrouping(setup, input_file_json)
    dataset = cal_grouping.get_dataset("calgroup001a", "default")
    assert dataset["p01"] == ["r001", "r002", "r003"]


def test_get_filelists(input_file_json):
    cal_grouping = CalGrouping(setup, input_file_json)
    files = cal_grouping.get_filelists("calgroup001a", "default", "dsp")
    assert len(files) == 3
    assert all(isinstance(file, Path) for file in files)
    assert files == [
        Path(setup["paths"]["tmp_filelists"]) / "all-l200-p01-r001-cal-dsp.filelist",
        Path(setup["paths"]["tmp_filelists"]) / "all-l200-p01-r002-cal-dsp.filelist",
        Path(setup["paths"]["tmp_filelists"]) / "all-l200-p01-r003-cal-dsp.filelist",
    ]


def test_get_par_files(input_file_json):
    cal_grouping = CalGrouping(setup, input_file_json)
    par_files = cal_grouping.get_par_files(catalog, "calgroup001a", "default", "dsp")
    assert isinstance(par_files, list)
    assert par_files == [
        str(
            Path(setup["paths"]["tmp_par"])
            / "l200-p01-r001-cal-20210101T000000Z-{channel}-par_dsp.yaml"
        ),
        str(
            Path(setup["paths"]["tmp_par"])
            / "l200-p01-r002-cal-20210102T000000Z-{channel}-par_dsp.yaml"
        ),
    ]


def test_get_plt_files(input_file_json):
    cal_grouping = CalGrouping(setup, input_file_json)
    plt_files = cal_grouping.get_plt_files(catalog, "calgroup001a", "default", "dsp")
    assert isinstance(plt_files, list)
    assert plt_files == [
        str(
            Path(setup["paths"]["tmp_plt"])
            / "l200-p01-r001-cal-20210101T000000Z-{channel}-plt_dsp.pkl"
        ),
        str(
            Path(setup["paths"]["tmp_plt"])
            / "l200-p01-r002-cal-20210102T000000Z-{channel}-plt_dsp.pkl"
        ),
    ]


def test_get_log_file(input_file_json):
    cal_grouping = CalGrouping(setup, input_file_json)
    log_file = cal_grouping.get_log_file(
        catalog, "calgroup001a", "default", "dsp", "timestamp", "par_dsp"
    )
    assert isinstance(log_file, str)


def test_get_timestamp(input_file_json):
    cal_grouping = CalGrouping(setup, input_file_json)
    timestamp = cal_grouping.get_timestamp(catalog, "calgroup001a", "default", "dsp")
    assert isinstance(timestamp, str)


def test_get_wildcard_constraints(input_file_json):
    cal_grouping = CalGrouping(setup, input_file_json)
    constraints = cal_grouping.get_wildcard_constraints("calgroup001a", "default")
    assert isinstance(constraints, str)


def test_get_dataset_unknown_partition(input_file_json):
    cal_grouping = CalGrouping(setup, input_file_json)
    with pytest.raises(KeyError, match="not defined for channel"):
        cal_grouping.get_dataset("no_such_partition", "default")


def test_get_dataset_missing_default(tmp_path):
    file = tmp_path / "input.json"
    file.write_text(json.dumps({"V01234A": {"part1": {"p01": ["r001"]}}}))
    cal_grouping = CalGrouping(setup, file)
    with pytest.raises(KeyError, match="no 'default' section"):
        cal_grouping.get_dataset("part1", "V01234A")


def test_expand_runs_malformed_range_yields_empty(tmp_path):
    # pylegendmeta's expand_runs returns [] for malformed ranges (no error raised)
    file = tmp_path / "input.json"
    file.write_text(json.dumps({"default": {"part": {"p01": "001..r005"}}}))
    cal_grouping = CalGrouping(setup, file)
    assert cal_grouping.datasets["default"]["part"]["p01"] == []

    file.write_text(json.dumps({"default": {"part": {"p01": ["r001..x005"]}}}))
    cal_grouping = CalGrouping(setup, file)
    assert cal_grouping.datasets["default"]["part"]["p01"] == []


def test_expand_runs_all_preserved(tmp_path):
    # "all" must not be expanded to []
    file = tmp_path / "input.json"
    file.write_text(json.dumps({"default": {"part": {"p01": "all"}}}))
    cal_grouping = CalGrouping(setup, file)
    assert cal_grouping.datasets["default"]["part"]["p01"] == "all"


def test_empty_partition_warns(input_file_json, caplog):
    cal_grouping = CalGrouping(setup, input_file_json)
    empty_catalog = Catalog.get([{"valid_from": "20210101T000000Z", "apply": []}])

    with caplog.at_level(logging.WARNING, logger="legenddataflow.methods.cal_grouping"):
        timestamp = cal_grouping.get_timestamp(
            empty_catalog, "calgroup001a", "default", "dsp"
        )
        log_file = cal_grouping.get_log_file(
            empty_catalog, "calgroup001a", "default", "dsp", "timestamp", "par_dsp"
        )

    assert timestamp == "20000101T000000Z"
    assert log_file == "log.log"
    assert caplog.text.count("resolved to no par files") == 2
