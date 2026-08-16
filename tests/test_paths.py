from __future__ import annotations

import logging
from pathlib import Path

import pytest
import yaml
from legenddataflowscripts.workflow import subst_vars

from legenddataflow.methods import paths

testprod = Path(__file__).parent / "dummy_cycle"

with (testprod / "config.yaml").open() as r:
    setup = yaml.safe_load(r)
subst_vars(setup, var_values={"_": str(testprod)})


def test_simple_resolvers():
    assert paths.tier_path(setup) == f"{testprod}/generated/tier"
    assert paths.pars_path(setup) == f"{testprod}/generated/par"
    assert paths.par_overwrite_path(setup) == f"{testprod}/inputs/dataprod/overrides"
    assert paths.metadata_path(setup) == f"{testprod}/inputs"
    assert paths.filelist_path(setup) == f"{testprod}/generated/tmp/filelists"


def test_sandbox_path():
    assert paths.sandbox_path(setup) == ""
    assert paths.sandbox_path({"paths": {}}) is None


def test_missing_path_key():
    with pytest.raises(KeyError, match=r"dataflow-config paths\.tier_daq is not set"):
        paths.tier_daq_path({"paths": {}})
    with pytest.raises(KeyError, match=r"dataflow-config paths\.tier_dsp is not set"):
        paths.get_tier_path({"paths": {}}, "dsp")


def test_get_tier_path():
    assert paths.get_tier_path(setup, "dsp") == f"{testprod}/generated/tier/dsp"
    assert paths.get_tier_path(setup, "skm") == f"{testprod}/generated/tier/skm"
    with pytest.raises(ValueError, match="no tier matching"):
        paths.get_tier_path(setup, "daq")
    with pytest.raises(ValueError, match="no tier matching"):
        paths.get_tier_path(setup, "not_a_tier")


def test_get_pars_path():
    assert paths.get_pars_path(setup, "hit") == f"{testprod}/generated/par/hit"
    with pytest.raises(ValueError, match="no tier matching"):
        paths.get_pars_path(setup, "skm")


def test_link_external_paths_creates_symlink(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    external = tmp_path / "other-prod" / "generated" / "tier" / "dsp"
    external.mkdir(parents=True)

    paths.link_external_paths({"paths": {"tier_dsp": str(external)}})

    default = tmp_path / "generated" / "tier" / "dsp"
    assert default.is_symlink()
    assert default.resolve() == external.resolve()

    # a second call with the same config is a no-op
    paths.link_external_paths({"paths": {"tier_dsp": str(external)}})
    assert default.resolve() == external.resolve()


def test_link_external_paths_removes_stale_symlink(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    default = tmp_path / "generated" / "tier" / "dsp"
    default.parent.mkdir(parents=True)
    default.symlink_to(tmp_path / "somewhere-else")

    paths.link_external_paths({"paths": {"tier_dsp": str(default)}})

    assert not default.is_symlink()
    assert not default.exists()


def test_link_external_paths_keeps_real_directory(tmp_path, monkeypatch, caplog):
    monkeypatch.chdir(tmp_path)
    default = tmp_path / "generated" / "tier" / "dsp"
    default.mkdir(parents=True)
    external = tmp_path / "other-prod" / "dsp"
    external.mkdir(parents=True)

    with caplog.at_level(logging.WARNING):
        paths.link_external_paths({"paths": {"tier_dsp": str(external)}})

    assert not default.is_symlink()
    assert "ignoring override" in caplog.text
