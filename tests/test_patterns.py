from __future__ import annotations

import copy
from pathlib import Path

import pytest
import yaml
from legenddataflowscripts.workflow import subst_vars

from legenddataflow.methods import patterns

testprod = Path(__file__).parent / "dummy_cycle"

with (testprod / "config.yaml").open() as r:
    setup = yaml.safe_load(r)
subst_vars(setup, var_values={"_": str(testprod)})


def test_key_patterns():
    assert (
        patterns.key_pattern() == "{experiment}-{period}-{run}-{datatype}-{timestamp}"
    )
    assert patterns.processing_pattern().endswith("-{processing_step}.{ext}")
    assert patterns.par_validity_pattern().startswith("{datatype}/{period}/{run}/")
    assert "{channel}" in patterns.full_channel_pattern()


def test_get_pattern_tier():
    assert (
        str(patterns.get_pattern_tier(setup, "dsp"))
        == f"{testprod}/generated/tier/dsp/{{datatype}}/{{period}}/{{run}}/"
        "{experiment}-{period}-{run}-{datatype}-{timestamp}-tier_dsp.lh5"
    )
    # concat patterns have no timestamp and no period/run directories
    assert (
        str(patterns.get_pattern_tier(setup, "evt_concat"))
        == f"{testprod}/generated/tier/evt/{{datatype}}/"
        "{experiment}-{period}-{run}-{datatype}-tier_evt.lh5"
    )
    # skm is phy-only and not per-run
    assert (
        str(patterns.get_pattern_tier(setup, "skm"))
        == f"{testprod}/generated/tier/skm/phy/"
        "{experiment}-{period}-{run}-{datatype}-tier_skm.lh5"
    )
    with pytest.raises(Exception, match="invalid tier"):
        patterns.get_pattern_tier(setup, "not_a_tier")


def test_get_pattern_tier_out_of_cycle():
    # the dummy config points tier_raw outside the production cycle
    out = patterns.get_pattern_tier(setup, "raw")
    assert str(out).startswith("/tmp/")
    # with the check disabled the configured path is returned unchanged
    kept = patterns.get_pattern_tier(setup, "raw", check_in_cycle=False)
    assert str(kept).startswith("/data2/")


def test_get_pattern_pars():
    assert (
        str(patterns.get_pattern_pars(setup, "dsp"))
        == f"{testprod}/generated/par/dsp/cal/{{period}}/{{run}}/"
        "{experiment}-{period}-{run}-cal-{timestamp}-par_dsp.yaml"
    )
    assert str(patterns.get_pattern_pars(setup, "dsp", name="eopt")).endswith(
        "par_dsp_eopt.yaml"
    )
    assert (
        str(patterns.get_pattern_pars(setup, "dsp", datatype=None)).count("{datatype}")
        == 2
    )
    with pytest.raises(Exception, match="invalid tier"):
        patterns.get_pattern_pars(setup, "skm")


def test_get_pattern_pars_out_of_cycle():
    modified = copy.deepcopy(setup)
    modified["paths"]["par_dsp"] = "/somewhere/else/par/dsp"
    out = patterns.get_pattern_pars(modified, "dsp", name="eopt", extension="pkl")
    assert str(out) == (
        "/tmp/{experiment}-{period}-{run}-cal-{timestamp}-par_dsp_eopt.pkl"
    )


def test_get_pattern_pars_tmp_channel():
    assert (
        str(patterns.get_pattern_pars_tmp_channel(setup, "dsp"))
        == f"{testprod}/generated/tmp/par/"
        "{experiment}-{period}-{run}-cal-{timestamp}-{channel}-par_dsp.yaml"
    )
    assert str(
        patterns.get_pattern_pars_tmp_channel(setup, "dsp", name="eopt")
    ).endswith("-{channel}-par_dsp_eopt.yaml")


def test_get_pattern_log():
    assert (
        str(patterns.get_pattern_log(setup, "build_dsp", "20230101T000000Z"))
        == f"{testprod}/generated/tmp/log/20230101T000000Z/build_dsp/"
        "{experiment}-{period}-{run}-{datatype}-{timestamp}-build_dsp.log"
    )
    assert str(
        patterns.get_pattern_log_channel(setup, "par_dsp", "20230101T000000Z")
    ).endswith("-cal-{timestamp}-{channel}-par_dsp.log")
