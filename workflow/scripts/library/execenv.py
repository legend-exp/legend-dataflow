# ruff: noqa: T201
from __future__ import annotations

import argparse
import os
import shutil
import string
import subprocess
from pathlib import Path

import yaml


def dataprod() -> None:
    """dataprod's command-line interface for installing and loading the software in the data production environment.

    .. code-block:: console

      $ dataprod --help
      $ dataprod load --help  # help section for a specific sub-command
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

    parser_load = subparsers.add_parser(
        "load", help="load data production environment and execute a given command"
    )
    parser_load.add_argument("config_file", help="production cycle configuration file", type=str)
    parser_load.add_argument(
        "command", help="command to run within the container", type=str, nargs="+"
    )
    parser_load.set_defaults(func=load)

    args = parser.parse_args()
    args.func(args)


def install(args) -> None:
    """
    This function installs user software in the data production environment.
    The software packages should be specified in the config.yaml file with the format:

    ```yaml
    setups:
        l200:
            pkg_versions:
                package_name: package_version
    ```
    """
    print(args.config_file)
    if not Path(args.config_file).is_file():
        msg = "config file is not a regular file"
        raise RuntimeError(msg)

    config_file_dir = Path(args.config_file).resolve().parent
    with Path(args.config_file).open() as r:
        config_dic = yaml.safe_load(r)

    exec_cmd = config_dic["setups"]["l200"]["execenv"]["cmd"]
    exec_arg = config_dic["setups"]["l200"]["execenv"]["arg"]
    path_src = config_dic["setups"]["l200"]["paths"]["src"]
    path_install = config_dic["setups"]["l200"]["paths"]["install"]
    path_cache = config_dic["setups"]["l200"]["paths"]["cache"]

    exec_cmd = string.Template(exec_cmd).substitute({"_": config_file_dir})
    exec_arg = string.Template(exec_arg).substitute({"_": config_file_dir})
    path_src = Path(string.Template(path_src).substitute({"_": config_file_dir}))
    path_install = Path(string.Template(path_install).substitute({"_": config_file_dir}))
    path_cache = Path(string.Template(path_cache).substitute({"_": config_file_dir}))

    if args.r:
        shutil.rmtree(path_install)
        shutil.rmtree(path_cache)

    pkg_list = ""
    for pkg, pkg_version in config_dic["setups"]["l200"]["pkg_versions"].items():
        if (path_src / pkg).exists():
            pkg_list += f" '{path_src / pkg}'"
        else:
            pkg_list += f" '{pkg_version}'"

    cmd_expr = (
        f"PYTHONUSERBASE={path_install} PIP_CACHE_DIR={path_cache} "
        f"{exec_cmd} {exec_arg} python3 -B -m pip install --no-warn-script-location {pkg_list}"
    )
    print("INFO: running:", cmd_expr)
    os.system(cmd_expr)


def load(args) -> None:
    """
    This function loads the data production environment and executes a given command.
    """

    if not Path(args.config_file).is_file():
        print("Error: config file does not exist")
        exit()

    config_file_dir = Path(args.config_file).resolve().parent
    with Path(args.config_file).open() as r:
        config_dic = yaml.safe_load(r)

    exec_cmd = config_dic["setups"]["l200"]["execenv"]["cmd"]
    exec_arg = config_dic["setups"]["l200"]["execenv"]["arg"]
    env_vars = config_dic["setups"]["l200"]["execenv"]["env"]
    path_install = config_dic["setups"]["l200"]["paths"]["install"]

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
