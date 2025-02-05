import os

import pytest
from dbetto import AttrsDict
from legenddataflow import execenv

os.environ["XDG_RUNTIME_DIR"] = "whatever"


@pytest.fixture(scope="module")
def config():
    return AttrsDict(
        {
            "paths": {"install": ".snakemake/software"},
            "execenv": {
                "cmd": "apptainer exec",
                "arg": "image.sif",
                "env": {
                    "VAR1": "val1",
                    "VAR2": "val2",
                },
            },
        }
    )


def test_execenv2str():
    assert (
        execenv._execenv2str(["cmd", "-v", "opt"], {"VAR1": "val1", "VAR2": "val2"})
        == "VAR1=val1 VAR2=val2 cmd -v opt"
    )


def test_execenv_prefix(config):
    cmd_expr, cmd_env = execenv.execenv_prefix(config, as_string=False)

    assert cmd_expr == [
        "apptainer",
        "exec",
        "--env=VAR1=val1",
        "--env=VAR2=val2",
        "--bind=whatever",
        "image.sif",
    ]
    assert cmd_env == config.execenv.env

    config.execenv.cmd = "docker run"
    cmd_expr, cmd_env = execenv.execenv_prefix(config, as_string=False)

    assert cmd_expr == [
        "docker",
        "run",
        "--env=VAR1=val1",
        "--env=VAR2=val2",
        "--volume=whatever:whatever",
        "image.sif",
    ]
    assert cmd_env == config.execenv.env

    config.execenv.cmd = "shifter"
    config.execenv.arg = "--image=legendexp/legend-base:latest"
    cmd_expr, cmd_env = execenv.execenv_prefix(config, as_string=False)

    assert cmd_expr == [
        "shifter",
        "--env=VAR1=val1",
        "--env=VAR2=val2",
        "--volume=whatever:whatever",
        "--image=legendexp/legend-base:latest",
    ]
    assert cmd_env == config.execenv.env

    cmd_str = execenv.execenv_prefix(config, as_string=True)
    assert cmd_str == (
        "VAR1=val1 VAR2=val2 "
        "shifter --env=VAR1=val1 --env=VAR2=val2 "
        "--volume=whatever:whatever "
        "--image=legendexp/legend-base:latest "
    )


def test_execenv_pyexe(config):
    cmd_str = execenv.execenv_pyexe(config, "dio-boe")

    assert cmd_str == (
        "VAR1=val1 VAR2=val2 "
        "shifter --env=VAR1=val1 --env=VAR2=val2 "
        "--volume=whatever:whatever "
        "--image=legendexp/legend-base:latest "
        ".snakemake/software/bin/dio-boe "
    )
