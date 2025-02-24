from __future__ import annotations

import argparse
import logging
import os
import shlex
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Iterable, Mapping

import colorlog
import dbetto
from dbetto import AttrsDict

from . import utils

log = logging.getLogger(__name__)


def _execenv2str(cmd_expr: Iterable, cmd_env: Mapping) -> str:
    return " ".join([f"{k}={v}" for k, v in cmd_env.items()]) + " " + " ".join(cmd_expr)


def apptainer_env_vars(cmdenv: Mapping) -> list[str]:
    return [f"--env={var}={val}" for var, val in cmdenv.items()]


def docker_env_vars(cmdenv: Mapping) -> list[str]:
    # same syntax
    return apptainer_env_vars(cmdenv)


def shifter_env_vars(cmdenv: Mapping) -> list[str]:
    # same syntax
    return apptainer_env_vars(cmdenv)


def execenv_prefix(
    config: AttrsDict, as_string: bool = True
) -> str | tuple[list, dict]:
    """Returns the software environment command prefix.

    For example: `apptainer run image.sif`

    Note
    ----
    If `as_string` is True, a space is appended to the returned string.
    """
    config = AttrsDict(config)

    cmdline = []
    cmdenv = {}
    if "execenv" in config and "env" in config.execenv:
        cmdenv |= config.execenv.env

    if "execenv" in config and "cmd" in config.execenv and "arg" in config.execenv:
        cmdline = shlex.split(config.execenv.cmd)

        has_xdg = False
        xdg_runtime_dir = os.getenv("XDG_RUNTIME_DIR")
        if xdg_runtime_dir:
            has_xdg = True

        if "env" in config.execenv:
            if any(exe in config.execenv.cmd for exe in ("apptainer", "singularity")):
                cmdline += apptainer_env_vars(config.execenv.env)
                if has_xdg:
                    cmdline += [f"--bind={xdg_runtime_dir}"]

            elif "docker" in config.execenv.cmd:
                cmdline += docker_env_vars(config.execenv.env)

            elif "shifter" in config.execenv.cmd:
                cmdline += shifter_env_vars(config.execenv.env)

            if (
                any(exe in config.execenv.cmd for exe in ("docker", "shifter"))
                and has_xdg
            ):
                cmdline += [f"--volume={xdg_runtime_dir}:{xdg_runtime_dir}"]

        # now we can add the arguments
        cmdline += shlex.split(config.execenv.arg)

    if as_string:
        return _execenv2str(cmdline, cmdenv) + " "

    return cmdline, cmdenv


def execenv_pyexe(
    config: AttrsDict, exename: str, as_string: bool = True
) -> str | tuple[list, dict]:
    """Returns the path to an executable installed in the virtualenv.

    For example: `apptainer run image.sif path/to/venv/bin/{exename}`

    Note
    ----
    If `as_string` is True, a space is appended to the returned string.
    """
    config = AttrsDict(config)

    cmdline, cmdenv = execenv_prefix(config, as_string=False)
    cmdline.append(f"{config.paths.install}/bin/{exename}")

    if as_string:
        return _execenv2str(cmdline, cmdenv) + " "

    return cmdline, cmdenv


def dataprod() -> None:
    """dataprod's CLI for installing and loading the software in the data production environment.

    .. code-block:: console

      $ dataprod --help
      $ dataprod install --help  # help section for a specific sub-command
    """

    parser = argparse.ArgumentParser(
        prog="dataprod", description="dataprod's command-line interface"
    )

    parser.add_argument(
        "-v", "--verbose", help="increase verbosity", action="store_true"
    )

    subparsers = parser.add_subparsers()
    parser_install = subparsers.add_parser(
        "install", help="install user software in data production environment"
    )
    parser_install.add_argument(
        "config_file", help="production cycle configuration file"
    )
    parser_install.add_argument(
        "-s",
        "--system",
        help="system running on",
        default="bare",
        type=str,
        required=False,
    )
    parser_install.add_argument(
        "-r",
        "--remove",
        help="remove software directory before installing software",
        action="store_true",
    )
    parser_install.add_argument(
        "-e",
        "--editable",
        help="install software with pip's --editable flag",
        action="store_true",
    )
    parser_install.set_defaults(func=install)

    parser_exec = subparsers.add_parser(
        "exec", help="load data production environment and execute a given command"
    )
    parser_exec.add_argument(
        "config_file", help="production cycle configuration file", type=str
    )
    parser_exec.add_argument(
        "--system", help="system running on", default="local", type=str, required=False
    )
    parser_exec.add_argument(
        "command", help="command to run within the container", type=str, nargs="+"
    )
    parser_exec.set_defaults(func=cmdexec)

    args = parser.parse_args()

    if args.verbose:
        handler = colorlog.StreamHandler()
        handler.setFormatter(
            colorlog.ColoredFormatter(
                "%(log_color)s%(name)s [%(levelname)s] %(message)s"
            )
        )

        logger = logging.getLogger("legenddataflow")
        logger.setLevel(logging.DEBUG)
        logger.addHandler(handler)

    args.func(args)


def install(args) -> None:
    """Installs user software in the data production environment.

    The software packages should be specified in the `config_file` with the
    format:

    ```yaml
    pkg_versions:
      - python_package_spec
      - ...
    ```

    .. code-block:: console

      $ dataprod install config.yaml
      $ dataprod install --editable config.yaml  # install legend-dataflow in editable mode
      $ dataprod install --remove config.yaml  # remove install directory
    """
    config_dict = AttrsDict(dbetto.utils.load_dict(args.config_file))
    config_loc = Path(args.config_file).resolve().parent

    utils.subst_vars(
        config_dict, var_values={"_": config_loc}, use_env=True, ignore_missing=False
    )
    config_dict["execenv"] = config_dict["execenv"][args.system]

    # path to virtualenv location
    path_install = config_dict.paths.install

    if args.remove and Path(path_install).exists():
        msg = f"removing: {path_install}"
        log.info(msg)
        shutil.rmtree(path_install)

    def _runcmd(cmd_expr, cmd_env, **kwargs):
        msg = "running: " + _execenv2str(cmd_expr, cmd_env)
        log.debug(msg)

        subprocess.run(cmd_expr, env=os.environ | cmd_env, check=True, **kwargs)

    cmd_prefix, cmd_env = execenv_prefix(config_dict, as_string=False)
    # HACK: get the full path to this python interpreter in case there is no execenv prefix
    python = sys.executable if cmd_prefix == [] else "python"
    python_venv, _ = execenv_pyexe(config_dict, "python", as_string=False)

    # we'll use uv from the virtualenv (installed below)
    uv_expr = [*python_venv, "-m", "uv", "--quiet"]

    # otherwise use python-venv
    cmd_expr = [*cmd_prefix, python, "-m", "venv", path_install]

    log.info(f"configuring virtual environment in {path_install}")
    _runcmd(cmd_expr, cmd_env)

    cmd_expr = [
        *python_venv,
        "-m",
        "pip",
        "--quiet",
        "--no-cache-dir",
        "install",
        "--upgrade",
        "--",
        "pip",
    ]

    log.info("upgrading pip")
    _runcmd(cmd_expr, cmd_env)

    # install uv
    cmd_expr = [
        *python_venv,
        "-m",
        "pip",
        "--quiet",
        "--no-cache-dir",
        "install",
        "--no-warn-script-location",
        "--",
        "uv",
    ]

    log.info("installing uv")
    _runcmd(cmd_expr, cmd_env)

    # and finally install legenddataflow with all dependencies
    # this must be done within the execenv, since jobs will be run within it

    cmd_expr = [
        *uv_expr,
        "pip",
        # "--no-cache",
        "install",
        "--prefix",
        path_install,
        str(config_loc),
    ]
    if args.editable:
        cmd_expr.insert(-1, "--editable")

    log.info("installing packages")
    _runcmd(cmd_expr, cmd_env)


def cmdexec(args) -> None:
    """Load the data production environment and execute a given command."""
    config_dict = AttrsDict(dbetto.utils.load_dict(args.config_file))
    config_loc = Path(args.config_file).resolve().parent

    utils.subst_vars(
        config_dict,
        var_values={"_": config_loc},
        use_env=True,
        ignore_missing=False,
    )
    config_dict["execenv"] = config_dict["execenv"][args.system]

    cmd_prefix, cmd_env = execenv_prefix(config_dict, as_string=False)
    cmd_expr = [*cmd_prefix, *args.command]

    msg = "running: " + _execenv2str(cmd_expr, cmd_env)
    log.debug(msg)

    subprocess.run(cmd_expr, env=os.environ | cmd_env, check=True)
