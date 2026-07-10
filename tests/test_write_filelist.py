from __future__ import annotations

import runpy
from types import SimpleNamespace

import snakemake.script

import legenddataflow.scripts.flow.write_filelist as write_filelist_module
from legenddataflow.scripts.flow.write_filelist import write_filelist


def test_write_filelist(tmp_path):
    out = tmp_path / "all-l200-p03-r000-phy-dsp.filelist"
    write_filelist(["/data/a.lh5", "/data/b.lh5"], out)
    assert out.read_text() == "/data/a.lh5\n/data/b.lh5\n"


def test_write_filelist_empty(tmp_path, capsys):
    out = tmp_path / "empty.filelist"
    write_filelist([], out, label="all-l200-p99-r999-phy")
    assert out.read_text() == ""
    captured = capsys.readouterr()
    assert "WARNING: no files found" in captured.out
    assert "all-l200-p99-r999-phy" in captured.out


def test_snakemake_glue(tmp_path, monkeypatch):
    # emulate a snakemake script: execution by injecting the snakemake object
    out = tmp_path / "out.filelist"
    fake = SimpleNamespace(
        input=["/data/a.lh5", "/data/b.lh5"],
        output=[str(out)],
        wildcards=SimpleNamespace(label="all-l200-p03-r000-phy"),
    )
    monkeypatch.setattr(snakemake.script, "snakemake", fake, raising=False)

    runpy.run_path(write_filelist_module.__file__)

    assert out.read_text() == "/data/a.lh5\n/data/b.lh5\n"
