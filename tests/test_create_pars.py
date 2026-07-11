from __future__ import annotations

from pathlib import Path

import pytest
import yaml
from dbetto import time
from legenddataflowscripts.workflow import subst_vars

from legenddataflow.methods import (
    FileKey,
    ParsKeyResolve,
    patterns,
)
from legenddataflow.methods.pars_loading import ParsCatalog

testprod = Path(__file__).parent / "dummy_cycle"

with (testprod / "config.yaml").open() as r:
    setup = yaml.safe_load(r)
subst_vars(setup, var_values={"_": str(testprod)})


def test_create_pars_keylist():
    key1 = FileKey("l200", "p00", "r000", "cal", "20230101T123456Z")
    assert ParsKeyResolve.entry_from_filekey(
        key1, {"cal": ["par_dsp"]}
    ).valid_from == time.unix_time("20230101T123456Z")
    key2 = FileKey("l200", "p00", "r000", "cal", "20230102T123456Z")
    assert ParsKeyResolve.match_keys(key1, key2) == key1
    key3 = FileKey("l200", "p00", "r000", "cal", "20230101T000000Z")
    assert ParsKeyResolve.match_keys(key1, key3) == key3
    assert ParsKeyResolve.generate_par_keylist([key1, key2, key3]) == [key3]
    pkey1 = ParsKeyResolve.entry_from_filekey(key1, {"cal": ["par_dsp"]})
    pkey2 = ParsKeyResolve.entry_from_filekey(
        FileKey("l200", "p00", "r000", "lar", "20230102T123456Z"), {"lar": "par_dsp"}
    )
    assert pkey2.file == [
        "lar/p00/r000/l200-p00-r000-lar-20230102T123456Z-par_dsp.yaml"
    ]
    ParsKeyResolve.match_entries(pkey1, pkey2)
    assert set(pkey2.file) == {
        "cal/p00/r000/l200-p00-r000-cal-20230101T123456Z-par_dsp.yaml",
        "lar/p00/r000/l200-p00-r000-lar-20230102T123456Z-par_dsp.yaml",
    }
    keylist = sorted(
        ParsKeyResolve.get_keys(
            "-*-*-*-cal",
            patterns.get_pattern_tier_daq(setup, extension="*", check_in_cycle=False),
        ),
        key=FileKey.get_unix_timestamp,
    )
    assert keylist == [
        FileKey("l200", "p00", "r000", "cal", "20230101T123456Z"),
        FileKey("l200", "p00", "r001", "cal", "20230202T004321Z"),
    ]

    keylist += ParsKeyResolve.get_keys(
        "-*-*-*-lar", patterns.get_pattern_tier_daq(setup, check_in_cycle=False)
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
        )[1].file
    ) == {
        "cal/p00/r000/l200-p00-r000-cal-20230101T123456Z-par_dsp.yaml",
        "lar/p00/r000/l200-p00-r000-lar-20230110T123456Z-par_dsp.yaml",
    }


def test_apply_run_override_missing_all(tmp_path):
    overwrite_validity = tmp_path / "validity.yaml"
    overwrite_validity.write_text(
        yaml.dump(
            [
                {
                    "valid_from": "20230101T000000Z",
                    "apply": ["a.yaml"],
                    "category": "other",
                }
            ]
        )
    )
    hit_catalog = ParsCatalog({"all": []})
    with pytest.raises(ValueError, match="has no 'all' system entry"):
        ParsKeyResolve.apply_run_override(hit_catalog, {}, overwrite_validity)

    overwrite_validity.write_text(
        yaml.dump([{"valid_from": "20230101T000000Z", "apply": ["a.yaml"]}])
    )
    hit_catalog = ParsCatalog({"other": []})
    with pytest.raises(ValueError, match="par catalog passed to apply_run_override"):
        ParsKeyResolve.apply_run_override(hit_catalog, {}, overwrite_validity)
