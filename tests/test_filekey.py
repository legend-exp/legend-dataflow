from __future__ import annotations

from pathlib import Path

import pytest
import yaml
from legenddataflowscripts.workflow import subst_vars

from legenddataflow.methods import ChannelProcKey, FileKey, paths, patterns

testprod = Path(__file__).parent / "dummy_cycle"

with (testprod / "config.yaml").open() as r:
    setup = yaml.safe_load(r)
subst_vars(setup, var_values={"_": str(testprod)})


def test_filekey():
    key = FileKey("l200", "p00", "r000", "cal", "20230101T123456Z")
    assert key.name == "l200-p00-r000-cal-20230101T123456Z"
    assert key._list() == ["l200", "p00", "r000", "cal", "20230101T123456Z"]
    keypart = "-l200-p00-r000-cal"
    key = FileKey.parse_keypart(keypart)
    assert key.name == "l200-p00-r000-cal-*"
    key = FileKey.from_string("l200-p00-r000-cal-20230101T123456Z")
    assert key.name == "l200-p00-r000-cal-20230101T123456Z"
    key = FileKey.get_filekey_from_filename(
        "l200-p00-r000-cal-20230101T123456Z-tier_dsp.lh5"
    )
    assert key.name == "l200-p00-r000-cal-20230101T123456Z"
    assert (
        key.get_path_from_filekey(patterns.get_pattern_tier(setup, "dsp"))[0]
        == f"{paths.get_tier_path(setup, 'dsp')}/cal/p00/r000/l200-p00-r000-cal-20230101T123456Z-tier_dsp.lh5"
    )
    assert (
        FileKey.get_filekey_from_pattern(
            key.get_path_from_filekey(patterns.get_pattern_tier(setup, "dsp"))[0],
            patterns.get_pattern_tier(setup, "dsp"),
        ).name
        == key.name
    )


def test_get_filekey_from_pattern_no_match():
    with pytest.raises(ValueError, match="does not match"):
        FileKey.get_filekey_from_pattern(
            "stray-file.txt", patterns.get_pattern_tier(setup, "dsp")
        )


def test_parse_keypart_invalid():
    with pytest.raises(ValueError, match="cannot be parsed as a FileKey keypart"):
        FileKey.parse_keypart("garbage")
    # keypart without the leading dash
    with pytest.raises(ValueError, match="cannot be parsed"):
        FileKey.parse_keypart("l200-p00-r000-cal")
    # ChannelProcKey keyparts must start with 'all'
    with pytest.raises(ValueError, match="ChannelProcKey keypart"):
        ChannelProcKey.parse_keypart("-l200-p00-r000-cal")


def test_get_path_from_filekey_dict_kwargs():
    key = FileKey("l200", "p00", "r000", "cal", "20230101T123456Z")
    # dict-valued kwargs are resolved via the key's field values
    assert key.get_path_from_filekey(
        "{experiment}-{run}-{extra}", extra={"cal": "x"}
    ) == ["l200-r000-x"]
    # non-intersecting dict kwargs are dropped instead of crashing
    assert key.get_path_from_filekey("{experiment}-{run}", extra={"nomatch": "x"}) == [
        "l200-r000"
    ]
