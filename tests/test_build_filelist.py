from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace

import pytest
import yaml

from legenddataflow.methods import FileKey, patterns
from legenddataflow.scripts.flow.build_filelist import (
    build_filelist,
    concat_phy_filenames,
    expand_runs,
    get_analysis_runs,
    get_filelist,
    get_keys,
    get_pattern,
)


@pytest.fixture
def config(tmp_path):
    paths = {"tier": str(tmp_path / "generated/tier")}
    for tier in ("daq", "raw", "dsp", "evt", "pet"):
        paths[f"tier_{tier}"] = str(tmp_path / f"generated/tier/{tier}")
    paths["tier_raw_blind"] = str(tmp_path / "generated/tier/raw-blind")
    paths["sandbox_path"] = str(tmp_path / "sandbox")
    return {"paths": paths}


def test_expand_runs():
    assert expand_runs({"phy": {"p01": "r001..r003"}}) == {
        "phy": {"p01": ["r001", "r002", "r003"]}
    }
    assert expand_runs({"phy": {"p01": ["r000", "r002..r004"]}}) == {
        "phy": {"p01": ["r000", "r002", "r003", "r004"]}
    }
    assert expand_runs({"cal": {"p02": "all"}}) == {"cal": {"p02": "all"}}


def test_get_analysis_runs_defaults():
    assert get_analysis_runs() == ({}, [])


def test_get_analysis_runs_yaml(tmp_path):
    runs_file = tmp_path / "runlists.yaml"
    runs_file.write_text(
        yaml.dump({"sel": {"phy": {"p03": "r000..r001"}}, "other": {}})
    )
    analysis_runs, ignore_keys = get_analysis_runs(
        analysis_runs_file=runs_file, file_selection="sel"
    )
    assert analysis_runs == {"phy": {"p03": ["r000", "r001"]}}
    assert ignore_keys == []

    with pytest.raises(ValueError, match="Unknown file selection"):
        get_analysis_runs(
            analysis_runs_file=runs_file, file_selection="not_a_selection"
        )


def test_get_analysis_runs_ignore_keylist(tmp_path):
    keylist = tmp_path / "ignore.keylist"
    keylist.write_text(
        "l200-p00-r000-cal-20230101T123456Z # some comment\n"
        "l200-p00-r001-cal-20230102T123456Z\n"
    )
    _, ignore_keys = get_analysis_runs(keylist)
    assert ignore_keys == [
        "l200-p00-r000-cal-20230101T123456Z",
        "l200-p00-r001-cal-20230102T123456Z",
    ]

    with pytest.raises(ValueError, match="no ignore_keys file found"):
        get_analysis_runs(tmp_path / "missing.yaml")

    bad = tmp_path / "ignore.txt"
    bad.write_text("key")
    with pytest.raises(ValueError, match="not in json, yaml or keylist format"):
        get_analysis_runs(bad)


def test_get_keys():
    keys = get_keys("-l200-p00-r000-cal-20230101T123456Z")
    assert len(keys) == 1
    assert keys[0].name == "l200-p00-r000-cal-20230101T123456Z"

    # `_`-separated multi-selectors expand to all combinations
    keys = get_keys("-l200-p00-r000_r001-cal_phy")
    assert sorted(k.name for k in keys) == [
        "l200-p00-r000-cal-*",
        "l200-p00-r000-phy-*",
        "l200-p00-r001-cal-*",
        "l200-p00-r001-phy-*",
    ]


def test_get_pattern(config):
    assert get_pattern(config, "dsp") == patterns.get_pattern_tier(
        config, "dsp", check_in_cycle=False
    )
    # tiers that search a different tier's files
    assert get_pattern(config, "blind") == patterns.get_pattern_tier(
        config, "raw", check_in_cycle=False
    )
    assert get_pattern(config, "skm") == patterns.get_pattern_tier(
        config, "pet", check_in_cycle=False
    )
    assert get_pattern(config, "evt_concat") == patterns.get_pattern_tier(
        config, "evt", check_in_cycle=False
    )


def test_concat_phy_filenames(config):
    pattern = patterns.get_pattern_tier(config, "evt", check_in_cycle=False)
    fnames = [
        str(pattern)
        .replace("{experiment}", "l200")
        .replace("{period}", "p03")
        .replace("{run}", run)
        .replace("{datatype}", "phy")
        .replace("{timestamp}", ts)
        for run, ts in [
            ("r000", "20230101T000000Z"),
            ("r000", "20230101T010000Z"),
            ("r001", "20230102T000000Z"),
        ]
    ]
    out = concat_phy_filenames(config, fnames, "evt_concat")
    # one output per run, with no timestamp component
    assert len(out) == 2
    assert out[0].endswith("l200-p03-r000-phy-tier_evt.lh5")
    assert out[1].endswith("l200-p03-r001-phy-tier_evt.lh5")


def _touch_raw_files(config, keys):
    raw_pattern = patterns.get_pattern_tier(config, "raw", check_in_cycle=False)
    for key in keys:
        f = Path(key.get_path_from_filekey(raw_pattern)[0])
        f.parent.mkdir(parents=True, exist_ok=True)
        f.touch()


def test_build_filelist(config):
    keys = [
        FileKey("l200", "p03", "r000", "cal", "20230101T000000Z"),
        FileKey("l200", "p03", "r000", "phy", "20230101T010000Z"),
        FileKey("l200", "p03", "r001", "phy", "20230102T000000Z"),
    ]
    _touch_raw_files(config, keys)
    search_pattern = patterns.get_pattern_tier(config, "raw", check_in_cycle=False)
    wildcard_keys = get_keys("-l200-p03-*-*")

    out = build_filelist(config, wildcard_keys, search_pattern, "dsp")
    assert len(out) == 3
    assert all("tier_dsp.lh5" in f for f in out)

    # ignore_keys drops the matching key
    out = build_filelist(
        config,
        wildcard_keys,
        search_pattern,
        "dsp",
        ignore_keys={"unprocessable": ["l200-p03-r000-cal-20230101T000000Z"]},
    )
    assert len(out) == 2
    assert all("phy" in f for f in out)

    # a flat list (as produced from .keylist files) works too
    out = build_filelist(
        config,
        wildcard_keys,
        search_pattern,
        "dsp",
        ignore_keys=["l200-p03-r000-cal-20230101T000000Z"],
    )
    assert len(out) == 2

    # analysis_runs restricts to the selected runs
    out = build_filelist(
        config,
        wildcard_keys,
        search_pattern,
        "dsp",
        analysis_runs={"phy": {"p03": ["r001"]}},
    )
    assert len(out) == 1
    assert "r001" in out[0]


def test_build_filelist_concat(config):
    keys = [
        FileKey("l200", "p03", "r000", "phy", "20230101T000000Z"),
        FileKey("l200", "p03", "r000", "phy", "20230101T010000Z"),
    ]
    evt_pattern = patterns.get_pattern_tier(config, "evt", check_in_cycle=False)
    for key in keys:
        f = Path(key.get_path_from_filekey(evt_pattern)[0])
        f.parent.mkdir(parents=True, exist_ok=True)
        f.touch()

    out = build_filelist(config, get_keys("-l200-p03-*-*"), evt_pattern, "evt_concat")
    # both timestamps of the run collapse into one concatenated target
    assert len(out) == 1
    assert out[0].endswith("l200-p03-r000-phy-tier_evt.lh5")


def test_get_filelist_wildcards(config):
    keys = [FileKey("l200", "p03", "r000", "cal", "20230101T000000Z")]
    _touch_raw_files(config, keys)
    search_pattern = patterns.get_pattern_tier(config, "raw", check_in_cycle=False)

    wildcards = SimpleNamespace(label="all-l200-p03-r000-cal", tier="dsp")
    out = get_filelist(wildcards, config, search_pattern)
    assert len(out) == 1
    assert out[0].endswith("l200-p03-r000-cal-20230101T000000Z-tier_dsp.lh5")
