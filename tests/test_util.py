import json
from pathlib import Path

from scripts.library import (
    CalibCatalog,
    FileKey,
    pars_catalog,
    pars_key_resolve,
    subst_vars,
    unix_time,
)
from scripts.library.patterns import get_pattern_tier_daq, get_pattern_tier_dsp
from scripts.library.utils import (
    par_dsp_path,
    par_overwrite_path,
    tier_dsp_path,
    tier_path,
)

testprod = Path(__file__).parent / "dummy_cycle"

with testprod.open() as r:
    setup = json.load(r)
subst_vars(setup, var_values={"_": str(testprod)})
setup = setup["setups"]["test"]


def test_util():
    assert tier_path(setup) == str(testprod / "generated/tier")
    assert unix_time("20230101T123456Z") == 1672572896.0


def test_filekey():
    key = FileKey("l200", "p00", "r000", "cal", "20230101T123456Z")
    assert key.name == "l200-p00-r000-cal-20230101T123456Z"
    assert key._list() == ["l200", "p00", "r000", "cal", "20230101T123456Z"]
    keypart = "-l200-p00-r000-cal"
    key = FileKey.parse_keypart(keypart)
    assert key.name == "l200-p00-r000-cal-*"
    key = FileKey.from_string("l200-p00-r000-cal-20230101T123456Z")
    assert key.name == "l200-p00-r000-cal-20230101T123456Z"
    key = FileKey.get_filekey_from_filename("l200-p00-r000-cal-20230101T123456Z-tier_dsp.lh5")
    assert key.name == "l200-p00-r000-cal-20230101T123456Z"
    assert (
        key.get_path_from_filekey(get_pattern_tier_dsp(setup))[0]
        == f"{tier_dsp_path(setup)}/cal/p00/r000/l200-p00-r000-cal-20230101T123456Z-tier_dsp.lh5"
    )
    assert (
        FileKey.get_filekey_from_pattern(
            key.get_path_from_filekey(get_pattern_tier_dsp(setup))[0],
            get_pattern_tier_dsp(setup),
        ).name
        == key.name
    )


def test_create_pars_keylist():
    key1 = FileKey("l200", "p00", "r000", "cal", "20230101T123456Z")
    assert (
        pars_key_resolve.from_filekey(key1, {"cal": ["par_dsp"]}).valid_from == "20230101T123456Z"
    )
    key2 = FileKey("l200", "p00", "r000", "cal", "20230102T123456Z")
    assert pars_key_resolve.match_keys(key1, key2) == key1
    key3 = FileKey("l200", "p00", "r000", "cal", "20230101T000000Z")
    assert pars_key_resolve.match_keys(key1, key3) == key3
    assert pars_key_resolve.generate_par_keylist([key1, key2, key3]) == [key3]
    pkey1 = pars_key_resolve.from_filekey(key1, {"cal": ["par_dsp"]})
    pkey2 = pars_key_resolve.from_filekey(
        FileKey("l200", "p00", "r000", "lar", "20230102T123456Z"), {"lar": ["par_dsp"]}
    )
    assert pkey2.apply == ["lar/p00/r000/l200-p00-r000-lar-20230102T123456Z-par_dsp.json"]
    pars_key_resolve.match_entries(pkey1, pkey2)
    assert set(pkey2.apply) == {
        "cal/p00/r000/l200-p00-r000-cal-20230101T123456Z-par_dsp.json",
        "lar/p00/r000/l200-p00-r000-lar-20230102T123456Z-par_dsp.json",
    }

    keylist = sorted(
        pars_key_resolve.get_keys("-*-*-*-cal", get_pattern_tier_daq(setup)),
        key=FileKey.get_unix_timestamp,
    )
    assert keylist == [
        FileKey("l200", "p00", "r000", "cal", "20230101T123456Z"),
        FileKey("l200", "p00", "r001", "cal", "20230202T004321Z"),
    ]

    keylist += pars_key_resolve.get_keys("-*-*-*-lar", get_pattern_tier_daq(setup))
    keylist = sorted(keylist, key=FileKey.get_unix_timestamp)
    assert keylist == [
        FileKey("l200", "p00", "r000", "cal", "20230101T123456Z"),
        FileKey("l200", "p00", "r000", "lar", "20230110T123456Z"),
        FileKey("l200", "p00", "r001", "cal", "20230202T004321Z"),
    ]

    pkeylist = pars_key_resolve.generate_par_keylist(keylist)
    assert pkeylist == keylist
    assert set(
        pars_key_resolve.match_all_entries(pkeylist, {"cal": ["par_dsp"], "lar": ["par_dsp"]})[
            1
        ].apply
    ) == {
        "cal/p00/r000/l200-p00-r000-cal-20230101T123456Z-par_dsp.json",
        "lar/p00/r000/l200-p00-r000-lar-20230110T123456Z-par_dsp.json",
    }


def test_pars_loading():
    pars_files = CalibCatalog.get_calib_files(
        Path(par_dsp_path(setup)) / "validity.jsonl", "20230101T123456Z"
    )
    assert pars_files == ["cal/p00/r000/l200-p00-r000-cal-20230101T123456Z-par_dsp.json"]

    par_override_files = CalibCatalog.get_calib_files(
        Path(par_overwrite_path(setup)) / "dsp" / "validity.jsonl", "20230101T123456Z"
    )

    pars_files, pars_files_overwrite = pars_catalog.match_pars_files(
        pars_files, par_override_files
    )

    assert pars_files == ["cal/p00/r000/l200-p00-r000-cal-20230101T123456Z-par_dsp.json"]

    assert set(pars_catalog.get_par_file(setup, "20230101T123456Z", "dsp")) == {
        (
            Path(par_dsp_path(setup))
            / "cal/p00/r000/l200-p00-r000-cal-20230101T123456Z-par_dsp.json",
        ),
        (
            Path(par_overwrite_path(setup))
            / "dsp/cal/p00/r000/l200-p00-r000-cal-T%-par_dsp_energy-overwrite.json",
        ),
    }
