import json
from pathlib import Path

from legenddataflow import (
    FileKey,
    ParsKeyResolve,
    patterns,
    subst_vars,
    utils,
)

testprod = Path(__file__).parent / "dummy_cycle"

with (testprod / "config.json").open() as r:
    setup = json.load(r)
subst_vars(setup, var_values={"_": str(testprod)})
setup = setup["setups"]["test"]


def test_util():
    assert utils.tier_path(setup) == str(testprod / "generated/tier")
    assert utils.unix_time("20230101T123456Z") == 1672572896.0


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
        == f"{utils.get_tier_path(setup, 'dsp')}/cal/p00/r000/l200-p00-r000-cal-20230101T123456Z-tier_dsp.lh5"
    )
    assert (
        FileKey.get_filekey_from_pattern(
            key.get_path_from_filekey(patterns.get_pattern_tier(setup, "dsp"))[0],
            utils.get_pattern_tier(setup, "dsp"),
        ).name
        == key.name
    )


def test_create_pars_keylist():
    key1 = FileKey("l200", "p00", "r000", "cal", "20230101T123456Z")
    assert (
        ParsKeyResolve.from_filekey(key1, {"cal": ["par_dsp"]}).valid_from
        == "20230101T123456Z"
    )
    key2 = FileKey("l200", "p00", "r000", "cal", "20230102T123456Z")
    assert ParsKeyResolve.match_keys(key1, key2) == key1
    key3 = FileKey("l200", "p00", "r000", "cal", "20230101T000000Z")
    assert ParsKeyResolve.match_keys(key1, key3) == key3
    assert ParsKeyResolve.generate_par_keylist([key1, key2, key3]) == [key3]
    pkey1 = ParsKeyResolve.from_filekey(key1, {"cal": ["par_dsp"]})
    pkey2 = ParsKeyResolve.from_filekey(
        FileKey("l200", "p00", "r000", "lar", "20230102T123456Z"), {"lar": "par_dsp"}
    )
    assert pkey2.apply == [
        "lar/p00/r000/l200-p00-r000-lar-20230102T123456Z-par_dsp.yaml"
    ]
    ParsKeyResolve.match_entries(pkey1, pkey2)
    assert set(pkey2.apply) == {
        "cal/p00/r000/l200-p00-r000-cal-20230101T123456Z-par_dsp.yaml",
        "lar/p00/r000/l200-p00-r000-lar-20230102T123456Z-par_dsp.yaml",
    }

    keylist = sorted(
        ParsKeyResolve.get_keys("-*-*-*-cal", patterns.get_pattern_tier_daq(setup)),
        key=FileKey.get_unix_timestamp,
    )
    assert keylist == [
        FileKey("l200", "p00", "r000", "cal", "20230101T123456Z"),
        FileKey("l200", "p00", "r001", "cal", "20230202T004321Z"),
    ]

    keylist += ParsKeyResolve.get_keys(
        "-*-*-*-lar", patterns.get_pattern_tier_daq(setup)
    )
    keylist = sorted(keylist, key=FileKey.get_unix_timestamp)
    assert keylist == [
        FileKey("l200", "p00", "r000", "cal", "20230101T123456Z"),
        FileKey("l200", "p00", "r000", "lar", "20230110T123456Z"),
        FileKey("l200", "p00", "r001", "cal", "20230202T004321Z"),
    ]

    pkeylist = ParsKeyResolve.generate_par_keylist(keylist)
    assert pkeylist == keylist
    assert set(
        ParsKeyResolve.match_all_entries(
            pkeylist, {"cal": ["par_dsp"], "lar": ["par_dsp"]}
        )[1].apply
    ) == {
        "cal/p00/r000/l200-p00-r000-cal-20230101T123456Z-par_dsp.json",
        "lar/p00/r000/l200-p00-r000-lar-20230110T123456Z-par_dsp.json",
    }
