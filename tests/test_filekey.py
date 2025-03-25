from pathlib import Path

import yaml
from legenddataflow import FileKey, paths, patterns, subst_vars

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
