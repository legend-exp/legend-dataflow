from __future__ import annotations

import pickle as pkl
import shelve

import pytest
import yaml

from legenddataflow.scripts.flow.merge_channels import merge_channels

CHANNELS = ["ch1000000", "ch1000001"]


def _channel_file(tmp_path, channel, ext):
    return tmp_path / f"l200-p03-r000-cal-20230101T000000Z-{channel}-par_dsp.{ext}"


def _run(monkeypatch, args):
    monkeypatch.setattr("sys.argv", ["merge-channels", *args])
    merge_channels()


def test_merge_channels_yaml(tmp_path, monkeypatch):
    infiles = []
    for i, channel in enumerate(CHANNELS):
        f = _channel_file(tmp_path, channel, "yaml")
        f.write_text(yaml.dump({"a": i}))
        infiles.append(str(f))
    out_file = tmp_path / "merged" / "l200-p03-r000-cal-20230101T000000Z-par_dsp.yaml"

    _run(monkeypatch, ["--input", *infiles, "--output", str(out_file)])

    out_dict = yaml.safe_load(out_file.read_text())
    assert out_dict == {"ch1000000": {"a": 0}, "ch1000001": {"a": 1}}


def test_merge_channels_extension_mismatch(tmp_path, monkeypatch):
    f = _channel_file(tmp_path, CHANNELS[0], "json")
    f.write_text("{}")
    out_file = tmp_path / "out.yaml"

    with pytest.raises(RuntimeError, match="extension does not match"):
        _run(monkeypatch, ["--input", str(f), "--output", str(out_file)])


def test_merge_channels_pickle(tmp_path, monkeypatch):
    infiles = []
    for i, channel in enumerate(CHANNELS):
        f = _channel_file(tmp_path, channel, "pkl")
        with f.open("wb") as w:
            pkl.dump({"result": i}, w)
        infiles.append(str(f))
    out_file = tmp_path / "l200-p03-r000-cal-20230101T000000Z-par_dsp.pkl"

    _run(monkeypatch, ["--input", *infiles, "--output", str(out_file)])

    with out_file.open("rb") as r:
        out_dict = pkl.load(r)
    assert out_dict == {"ch1000000": {"result": 0}, "ch1000001": {"result": 1}}


def test_merge_channels_shelve(tmp_path, monkeypatch):
    infiles = []
    for i, channel in enumerate(CHANNELS):
        f = _channel_file(tmp_path, channel, "pkl")
        with f.open("wb") as w:
            pkl.dump({"result": i, "common": {"shared": True}}, w)
        infiles.append(str(f))
    out_file = tmp_path / "l200-p03-r000-cal-20230101T000000Z-plt_dsp.dir"

    _run(monkeypatch, ["--input", *infiles, "--output", str(out_file)])

    with shelve.open(str(out_file.with_suffix(""))) as shelf:
        assert shelf["ch1000000"] == {"result": 0}
        assert shelf["ch1000001"] == {"result": 1}
        # per-channel "common" blocks are collected into one entry
        assert shelf["common"] == {
            "ch1000000": {"shared": True},
            "ch1000001": {"shared": True},
        }
