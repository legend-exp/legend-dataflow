from __future__ import annotations

import pickle as pkl

import lh5
import pytest
import yaml
from lgdo.types import Array, Struct

from legenddataflow.scripts.flow.merge_in_channel import merge_in_channel

CHANNEL = "ch1000000"


def _infile(tmp_path, step, ext):
    return tmp_path / f"l200-p03-r000-cal-20230101T000000Z-{CHANNEL}-{step}.{ext}"


@pytest.fixture
def run(monkeypatch):
    def _run(args):
        monkeypatch.setattr("sys.argv", ["merge-in-channels", *args])
        merge_in_channel()

    return _run


def test_merge_in_channel_yaml(tmp_path, run):
    f1 = _infile(tmp_path, "par_dsp_eopt", "yaml")
    f1.write_text(yaml.dump({"eopt": {"a": 1}}))
    f2 = _infile(tmp_path, "par_dsp_dplms", "yaml")
    f2.write_text(yaml.dump({"dplms": {"b": 2}}))
    out_file = tmp_path / "merged" / f"l200-p03-r000-cal-{CHANNEL}-par_dsp.yaml"

    run(["--input", str(f1), str(f2), "--output", str(out_file)])

    assert yaml.safe_load(out_file.read_text()) == {
        "eopt": {"a": 1},
        "dplms": {"b": 2},
    }


def test_merge_in_channel_pickle(tmp_path, run):
    f1 = _infile(tmp_path, "par_dsp_eopt", "pkl")
    with f1.open("wb") as w:
        pkl.dump({"eopt": 1}, w)
    f2 = _infile(tmp_path, "par_dsp_dplms", "pkl")
    with f2.open("wb") as w:
        pkl.dump({"dplms": 2}, w)
    out_file = tmp_path / f"l200-p03-r000-cal-{CHANNEL}-par_dsp.pkl"

    run(["--input", str(f1), str(f2), "--output", str(out_file)])

    with out_file.open("rb") as r:
        assert pkl.load(r) == {"eopt": 1, "dplms": 2}


def test_merge_in_channel_lh5(tmp_path, run):
    f1 = _infile(tmp_path, "par_dsp_eopt", "lh5")
    lh5.write(
        Struct({"eopt": Array([1, 2, 3])}), name=CHANNEL, lh5_file=str(f1), wo_mode="a"
    )
    f2 = _infile(tmp_path, "par_dsp_dplms", "lh5")
    lh5.write(
        Struct({"dplms": Array([4, 5])}), name=CHANNEL, lh5_file=str(f2), wo_mode="a"
    )
    out_file = tmp_path / f"l200-p03-r000-cal-20230101T000000Z-{CHANNEL}-par_dsp.lh5"

    run(["--input", str(f1), str(f2), "--output", str(out_file)])

    res = lh5.read(CHANNEL, str(out_file))
    # the identifier of each input's processing step becomes a struct field
    assert set(res.keys()) == {"eopt", "dplms"}
    assert list(res["eopt"]) == [1, 2, 3]
    assert list(res["dplms"]) == [4, 5]


def test_merge_in_channel_lh5_multipart_identifier(tmp_path, run):
    # identifiers may contain underscores and must be kept whole
    f = _infile(tmp_path, "par_dsp_svm_train", "lh5")
    lh5.write(
        Struct({"svm_train": Array([1, 2])}), name=CHANNEL, lh5_file=str(f), wo_mode="a"
    )
    out_file = tmp_path / f"l200-p03-r000-cal-20230101T000000Z-{CHANNEL}-par_dsp.lh5"

    run(["--input", str(f), "--output", str(out_file)])

    res = lh5.read(CHANNEL, str(out_file))
    assert set(res.keys()) == {"svm_train"}


def test_merge_in_channel_lh5_no_identifier(tmp_path, run):
    f = _infile(tmp_path, "par_dsp", "lh5")
    lh5.write(Struct({"a": Array([1])}), name=CHANNEL, lh5_file=str(f), wo_mode="a")
    out_file = tmp_path / f"l200-p03-r000-cal-20230101T000000Z-{CHANNEL}-par_dsp.lh5"

    with pytest.raises(ValueError, match="has no identifier part"):
        run(["--input", str(f), "--output", str(out_file)])


def test_merge_in_channel_lh5_mixed_channels(tmp_path, run):
    f1 = _infile(tmp_path, "par_dsp_eopt", "lh5")
    lh5.write(Struct({"eopt": Array([1])}), name=CHANNEL, lh5_file=str(f1), wo_mode="a")
    other = "ch1000001"
    f2 = tmp_path / f"l200-p03-r000-cal-20230101T000000Z-{other}-par_dsp_dplms.lh5"
    lh5.write(Struct({"dplms": Array([2])}), name=other, lh5_file=str(f2), wo_mode="a")
    out_file = tmp_path / f"l200-p03-r000-cal-20230101T000000Z-{CHANNEL}-par_dsp.lh5"

    with pytest.raises(ValueError, match="span multiple channels"):
        run(["--input", str(f1), str(f2), "--output", str(out_file)])


def test_merge_in_channel_lh5_missing_field(tmp_path, run):
    f = _infile(tmp_path, "par_dsp_eopt", "lh5")
    lh5.write(Struct({"wrong": Array([1])}), name=CHANNEL, lh5_file=str(f), wo_mode="a")
    out_file = tmp_path / f"l200-p03-r000-cal-20230101T000000Z-{CHANNEL}-par_dsp.lh5"

    with pytest.raises(KeyError, match="has no field 'eopt'"):
        run(["--input", str(f), "--output", str(out_file)])
