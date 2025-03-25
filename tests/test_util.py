import os
from pathlib import Path

import pytest
import yaml
from legenddataflow.paths import tier_path
from legenddataflow.utils import (
    subst_vars,
    subst_vars_impl,
    subst_vars_in_snakemake_config,
)

testprod = Path(__file__).parent / "dummy_cycle"
config_filename = testprod / "config.yaml"

with (config_filename).open() as r:
    setup = yaml.safe_load(r)
subst_vars(setup, var_values={"_": str(testprod)})


@pytest.fixture
def mock_os_environ(monkeypatch):
    monkeypatch.setenv("PRODENV", "prod")
    return os.environ


class mock_workflow_class:
    def __init__(self):
        self.overwrite_configfiles = [config_filename]


@pytest.fixture(scope="module")
def mock_workflow():
    return mock_workflow_class()


def test_util():
    assert tier_path(setup) == str(testprod / "generated/tier")


def test_subst_vars_impl():
    input_str = "$PRODENV"
    var_values = {"PRODENV": "prod", "_": "/path/to/prod"}
    # test string substitution
    result = subst_vars_impl(input_str, var_values)
    assert result == "prod"

    input_str = "$PROD"
    # test string substitution with ignore missing
    result = subst_vars_impl(input_str, var_values, ignore_missing=True)
    assert result == input_str
    with pytest.raises(KeyError):
        subst_vars_impl(input_str, var_values, ignore_missing=False)

    # test dictionary substitution
    input_dict = {"prodenv": "$PRODENV", "tier": "$_/tier"}
    expected = {"prodenv": "prod", "tier": "/path/to/prod/tier"}
    result = subst_vars_impl(input_dict, var_values)
    assert result == expected

    # test list substitution
    input_list = ["$PRODENV in $_"]
    expected = ["prod in /path/to/prod"]
    result = subst_vars_impl(input_list, var_values)
    assert result == expected

    # Test mixed structure with dictionary and list
    input_dict = {
        "prodenv": "$PRODENV",
        "path": {"tier": "$_/tier"},
        "test": ["$PRODENV in $_"],
    }
    expected = {
        "prodenv": "prod",
        "path": {"tier": "/path/to/prod/tier"},
        "test": ["prod in /path/to/prod"],
    }
    result = subst_vars_impl(input_dict, var_values)
    assert result == expected

    input_data = [42, 3.14, None]
    for inval in input_data:
        result = subst_vars_impl(inval, {})
        assert result == inval


def test_subst_vars(mock_os_environ):  # noqa: ARG001
    # test no var values provided and no env
    props = {
        "prodenv": "$PRODENV",
        "path": {"tier": "$_/tier"},
        "test": ["$PRODENV in $_"],
    }
    expected = {
        "prodenv": "$PRODENV",
        "path": {"tier": "$_/tier"},
        "test": ["$PRODENV in $_"],
    }
    with pytest.raises(KeyError):
        subst_vars(props, ignore_missing=False)
    subst_vars(props, ignore_missing=True)
    assert props == expected

    # test var values provided and no env
    props = {
        "prodenv": "$PRODENV",
        "path": {"tier": "$_/tier"},
        "test": ["$PRODENV in $_"],
    }
    var_values = {"PRODENV": "prod", "_": "/path/to/prod"}
    expected = {
        "prodenv": "prod",
        "path": {"tier": "/path/to/prod/tier"},
        "test": ["prod in /path/to/prod"],
    }
    subst_vars(props, var_values)
    subst_vars(props, var_values)
    assert props == expected

    # test no var values provided and env
    props = {
        "prodenv": "$PRODENV",
    }
    expected = {
        "prodenv": "prod",
    }
    subst_vars(props, use_env=True)
    assert props == expected

    # test var values provided and env
    props = {
        "prodenv": "$PRODENV",
        "path": {"tier": "$_/tier"},
        "test": ["$PRODENV in $_"],
    }
    var_values = {"_": "/path/to/prod"}
    expected = {
        "prodenv": "prod",
        "path": {"tier": "/path/to/prod/tier"},
        "test": ["prod in /path/to/prod"],
    }
    subst_vars(props, var_values, use_env=True)
    assert props == expected


def test_subst_vars_in_snakemake_config(mock_workflow, mock_os_environ):  # noqa: ARG001
    setup = {
        "prodenv": "$PRODENV",
        "path": {"tier": "$_/generated/tier"},
        "execenv": {
            "bare": {
                "env": {
                    "VAR1": "val1",
                }
            },
            "lngs": {
                "cmd": "apptainer exec",
                "arg": "$PRODENV/container.sif",
                "env": {
                    "VAR2": "val2",
                },
            },
        },
    }
    subst_vars_in_snakemake_config(mock_workflow, setup)
    assert setup["prodenv"] == "prod"
    assert setup["path"]["tier"] == str(testprod / "generated/tier")
    assert "lngs" not in setup["execenv"]
    assert setup["execenv"]["env"] == {"VAR1": "val1"}
    assert "cmd" not in setup["execenv"]
    assert "arg" not in setup["execenv"]
    setup = {
        "prodenv": "$PRODENV",
        "path": {"tier": "$_/generated/tier"},
        "execenv": {
            "bare": {
                "env": {
                    "VAR1": "val1",
                }
            },
            "lngs": {
                "cmd": "apptainer exec",
                "arg": "$PRODENV/container.sif",
                "env": {
                    "VAR2": "val2",
                },
            },
        },
        "system": "lngs",
    }
    subst_vars_in_snakemake_config(mock_workflow, setup)
    assert setup["prodenv"] == "prod"
    assert setup["path"]["tier"] == str(testprod / "generated/tier")
    assert "local" not in setup["execenv"]
    assert setup["execenv"]["env"] == {"VAR2": "val2"}
    assert setup["execenv"]["cmd"] == "apptainer exec"
    assert setup["execenv"]["arg"] == "prod/container.sif"
