from __future__ import annotations

import argparse
import logging
import os
import shlex
import shutil
import string
import subprocess
from pathlib import Path

import dbetto
from dbetto import AttrsDict
from packaging.requirements import Requirement

from . import utils

log = logging.getLogger(__name__)


def execenv_prefix(config, aslist=False):
    """Returns the software environment command prefix.

    For example: `apptainer run image.sif`
    """
    config = AttrsDict(config)

    cmdline = shlex.split(config.execenv.cmd)
    if "env" in config.execenv:
        cmdline += [f"--env={var}={val}" for var, val in config.execenv.env.items()]

    cmdline += shlex.split(config.execenv.arg)

    if aslist:
        return cmdline
    return " ".join(cmdline)


def execenv_python(config, aslist=False):
    """Returns the Python interpreter command.

    For example: `apptainer run image.sif python`
    """
    config = AttrsDict(config)

    cmdline = execenv_prefix(config, aslist=True)
    cmdline.append(f"{config.paths.install}/bin/python")

    if aslist:
        return cmdline
    return " ".join(cmdline)


def execenv_smk_py_script(workflow, config, scriptname, aslist=False):
    """Returns the command used to run a Python script for a Snakemake rule.

    For example: `apptainer run image.sif python path/to/script.py`
    """
    config = AttrsDict(config)

    cmdline = execenv_python(config, aslist=True)
    cmdline.append(f"{workflow.basedir}/scripts/{scriptname}")

    if aslist:
        return cmdline
    return " ".join(cmdline)


def dataprod() -> None:
    """dataprod's command-line interface for installing and loading the software in the data production environment.

    .. code-block:: console

      $ dataprod --help
      $ dataprod exec --help  # help section for a specific sub-command
    """

    parser = argparse.ArgumentParser(
        prog="dataprod", description="dataprod's command-line interface"
    )

    subparsers = parser.add_subparsers()
    parser_install = subparsers.add_parser(
        "install", help="install user software in data production environment"
    )
    parser_install.add_argument(
        "config_file", help="production cycle configuration file", type=str
    )
    parser_install.add_argument(
        "-r", help="remove software directory before installing software", action="store_true"
    )
    parser_install.set_defaults(func=install)

    parser_exec = subparsers.add_parser(
        "exec", help="load data production environment and execute a given command"
    )
    parser_exec.add_argument("config_file", help="production cycle configuration file", type=str)
    parser_exec.add_argument(
        "command", help="command to run within the container", type=str, nargs="+"
    )
    parser_exec.set_defaults(func=cmdexec)

    args = parser.parse_args()
    args.func(args)


def install(args) -> None:
    """
    This function installs user software in the data production environment.
    The software packages should be specified in the config.yaml file with the
    format:

    ```yaml
    pkg_versions:
      - python_package_spec
      - ...
    ```
    """
    config_dict = AttrsDict(dbetto.utils.load_dict(args.config_file))
    config_loc = Path(args.config_file).resolve().parent
    path_install = config_dict.paths.install

    if args.r and Path(path_install).exists():
        shutil.rmtree(path_install)

    utils.subst_vars(
        config_dict,
        var_values={"_": config_loc},
        use_env=True,
        ignore_missing=False,
    )

    cmd_env = {}

    def _runcmd(cmd_env, cmd_expr):
        msg = (
            "running:"
            + " ".join([f"{k}={v}" for k, v in cmd_env.items()])
            + " "
            + " ".join(cmd_expr),
        )
        log.debug(msg)

        subprocess.run(cmd_expr, env=cmd_env, check=True)

    # configure venv
    cmd_expr = [*execenv_prefix(config_dict, aslist=True), "python3", "-m", "venv", path_install]

    log.info(f"configuring virtual environment in {path_install}")
    _runcmd(cmd_env, cmd_expr)

    cmd_expr = [
        *execenv_python(config_dict, aslist=True),
        "-m",
        "pip",
        "--no-cache-dir",
        "install",
        "--upgrade",
        "pip",
    ]

    log.info("upgrading pip")
    _runcmd(cmd_env, cmd_expr)

    # install uv
    cmd_expr = [
        *execenv_python(config_dict, aslist=True),
        "-m",
        "pip",
        "--no-cache-dir",
        "install",
        "--no-warn-script-location",
        "uv",
    ]

    log.info("installing uv")
    _runcmd(cmd_env, cmd_expr)

    # now packages

    path_src = Path(config_dict.paths.src)
    pkg_list = []
    for spec in config_dict.pkg_versions:
        pkg = Requirement(spec).name
        if (path_src / pkg).exists():
            pkg_list.append(str(path_src / pkg))
        else:
            pkg_list.append(spec)

    cmd_base = [
        *execenv_python(config_dict, aslist=True),
        "-m",
        "uv",
        "pip",
        "--no-cache",
        "install",
    ]

    cmd_expr = cmd_base + pkg_list

    log.info("installing packages")
    _runcmd(cmd_env, cmd_expr)

    # and finally legenddataflow

    cmd_expr = [
        *execenv_python(config_dict, aslist=True),
        "-m",
        "uv",
        "pip",
        "--no-cache",
        "install",
        # "--editable",  # TODO do we really want this?
        str(config_loc),
    ]

    log.info("installing packages")
    _runcmd(cmd_env, cmd_expr)


def cmdexec(args) -> None:
    """
    This function loads the data production environment and executes a given command.
    """
    config_file_dir = Path(args.config_file).resolve().parent
    config_dict = AttrsDict(dbetto.utils.load_dict(args.config_file))

    exec_cmd = config_dict.execenv.cmd
    exec_arg = config_dict.execenv.arg
    env_vars = config_dict.execenv.env
    path_install = config_dict.paths.install

    exec_cmd = string.Template(exec_cmd).substitute({"_": config_file_dir})
    exec_arg = string.Template(exec_arg).substitute({"_": config_file_dir})
    path_install = string.Template(path_install).substitute({"_": config_file_dir})

    xdg_runtime_dir = os.getenv("XDG_RUNTIME_DIR")
    if xdg_runtime_dir:
        subprocess.run(
            [*(exec_cmd.split()), exec_arg, *args.command],
            env=dict(
                PYTHONUSERBASE=path_install,
                APPTAINERENV_APPEND_PATH=f":{path_install}/bin",
                APPTAINER_BINDPATH=xdg_runtime_dir,
                **env_vars,
            ),
            check=True,
        )
    else:
        subprocess.run(
            [*(exec_cmd.split()), exec_arg, *args.command],
            env=dict(
                PYTHONUSERBASE=path_install,
                APPTAINERENV_APPEND_PATH=f":{path_install}/bin",
                APPTAINER_BINDPATH=xdg_runtime_dir,
                **env_vars,
            ),
            check=True,
        )
