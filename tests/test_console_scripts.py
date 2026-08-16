"""Smoke tests for the console-script entry points in ``[project.scripts]``.

Each entry point is imported and invoked with ``--help``. This catches broken
imports (e.g. moved third-party symbols), misconfigured entry points and
argparse setup errors without needing any input data.
"""

from __future__ import annotations

import importlib
import sys
import tomllib
from pathlib import Path

import pytest

with (Path(__file__).parent.parent / "pyproject.toml").open("rb") as f:
    console_scripts = tomllib.load(f)["project"]["scripts"]


# the heavy scientific imports (pygama, matplotlib, ...) triggered by these
# modules may emit warnings; they are not what this smoke test is about
@pytest.mark.filterwarnings("ignore")
@pytest.mark.parametrize("script", sorted(console_scripts))
def test_console_script_help(script, monkeypatch, capsys):
    module_name, func_name = console_scripts[script].split(":")

    module = importlib.import_module(module_name)
    func = getattr(module, func_name)
    assert callable(func)

    monkeypatch.setattr(sys, "argv", [script, "--help"])
    with pytest.raises(SystemExit) as excinfo:
        func()
    assert excinfo.value.code == 0
    assert "usage" in capsys.readouterr().out.lower()
