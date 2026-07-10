from __future__ import annotations

import runpy
from types import SimpleNamespace

import pytest
import snakemake.script

import legenddataflow.scripts.flow.complete_run as complete_run_module
from legenddataflow.scripts.flow.complete_run import (
    build_file_db_config,
    check_log_files,
    fformat,
    find_gen_runs,
)


@pytest.fixture
def setup(tmp_path):
    paths = {"tier": str(tmp_path / "generated/tier")}
    for tier in ("raw", "dsp", "evt"):
        paths[f"tier_{tier}"] = str(tmp_path / f"generated/tier/{tier}")
    return {
        "paths": paths,
        "table_format": {"raw": "ch{ch:07d}/raw", "evt": "{grp}/evt"},
    }


def _make_logs(log_path):
    sub = log_path / "20230101T000000Z" / "build_dsp"
    sub.mkdir(parents=True)
    (sub / "bad.log").write_text(
        "INFO: starting\nERROR: something broke\nWARNING: be careful\n"
    )
    (sub / "clean.log").write_text("INFO: all fine\n")
    return sub


def test_check_log_files_with_warning_file(tmp_path):
    log_path = tmp_path / "log"
    _make_logs(log_path)
    summary = tmp_path / "summary" / "summary.log"
    warning = tmp_path / "summary" / "warning.log"

    check_log_files(log_path, summary, "all-l200.gen", warning_file=warning)

    summary_text = summary.read_text()
    assert "generated at" in summary_text
    assert "with errors" in summary_text
    assert "bad.log : ERROR: something broke" in summary_text
    warning_text = warning.read_text()
    assert "with warnings" in warning_text
    assert "bad.log : WARNING: be careful" in warning_text
    # log files are consumed and empty directories cleaned up
    assert not log_path.exists()


def test_check_log_files_clean(tmp_path):
    log_path = tmp_path / "log"
    sub = log_path / "build_dsp"
    sub.mkdir(parents=True)
    (sub / "clean.log").write_text("INFO: all fine\n")
    summary = tmp_path / "summary.log"
    warning = tmp_path / "warning.log"

    check_log_files(log_path, summary, "all-l200.gen", warning_file=warning)

    assert "with no errors" in summary.read_text()
    assert "with no warnings" in warning.read_text()
    assert not log_path.exists()


def test_check_log_files_no_warning_file(tmp_path):
    log_path = tmp_path / "log"
    _make_logs(log_path)
    summary = tmp_path / "summary.log"

    check_log_files(log_path, summary, "all-l200.gen")

    summary_text = summary.read_text()
    assert "with errors" in summary_text
    assert "ERROR: something broke" in summary_text
    assert not log_path.exists()


def test_find_gen_runs(tmp_path):
    # parent path contains a hyphen on purpose: run parsing for concat tiers
    # must use the filename, not the full path
    gen = tmp_path / "prod-blind" / "generated" / "tier"
    (gen / "dsp" / "phy" / "p03" / "r000").mkdir(parents=True)
    (gen / "evt" / "phy").mkdir(parents=True)
    (gen / "evt" / "phy" / "l200-p03-r001-phy-tier_evt.lh5").touch()

    assert find_gen_runs(gen) == {"phy/p03/r000", "phy/p03/r001"}


def test_fformat(setup):
    out = fformat(setup, "raw")
    # tier path prefix is stripped, leaving the relative file pattern
    assert out.startswith("/{datatype}/{period}/{run}/")
    assert out.endswith("-tier_raw.lh5")
    assert setup["paths"]["tier_raw"] not in out


def test_build_file_db_config_without_prodenv(setup, tmp_path, monkeypatch):
    monkeypatch.delenv("PRODENV", raising=False)

    cfg = build_file_db_config(setup, tmp_path / "filedb")

    assert cfg["data_dir"] == "/"
    assert cfg["tier_dirs"] == {
        "raw": setup["paths"]["tier_raw"],
        "evt": setup["paths"]["tier_evt"],
    }
    assert cfg["table_format"] == setup["table_format"]
    assert cfg["file_format"]["raw"].startswith("/{datatype}")


def test_build_file_db_config_with_prodenv(setup, tmp_path, monkeypatch):
    monkeypatch.setenv("PRODENV", str(tmp_path))

    cfg = build_file_db_config(setup, tmp_path / "filedb")

    assert cfg["data_dir"] == "$PRODENV"
    # tier dirs are relative to $PRODENV
    assert cfg["tier_dirs"]["raw"] == "/generated/tier/raw"
    assert cfg["tier_dirs"]["evt"] == "/generated/tier/evt"


def test_snakemake_glue(tmp_path, monkeypatch):
    # emulate a snakemake script: execution by injecting the snakemake object;
    # FileDB building is disabled so only the log check and touch run
    log_path = tmp_path / "log"
    _make_logs(log_path)
    setup = {
        "paths": {"tmp_log": str(log_path)},
        "build_file_dbs": False,
        "check_log_files": True,
    }
    fake = SimpleNamespace(
        params=SimpleNamespace(
            setup=setup, filedb_path=str(tmp_path / "filedb"), ignore_keys_file=None
        ),
        wildcards=SimpleNamespace(tier="dsp", label="all-l200-p03-r000-phy-dsp"),
        threads=1,
        output=SimpleNamespace(
            summary_log=str(tmp_path / "summary.log"),
            warning_log=str(tmp_path / "warning.log"),
            gen_output=str(tmp_path / "all-l200-p03-r000-phy-dsp.gen"),
        ),
    )
    monkeypatch.setattr(snakemake.script, "snakemake", fake, raising=False)

    runpy.run_path(complete_run_module.__file__)

    assert "with errors" in (tmp_path / "summary.log").read_text()
    assert "with warnings" in (tmp_path / "warning.log").read_text()
    assert (tmp_path / "all-l200-p03-r000-phy-dsp.gen").exists()
    assert not log_path.exists()
