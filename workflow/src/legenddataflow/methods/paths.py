"""
This module contains resolvers for the config.json dictionary
"""

from __future__ import annotations

import logging
import os
from pathlib import Path

import yaml
from legenddataflowscripts.workflow import subst_vars


def sandbox_path(setup):
    if "sandbox_path" in setup["paths"]:
        return setup["paths"]["sandbox_path"]
    return None


def tier_daq_path(setup):
    return setup["paths"]["tier_daq"]


def tier_raw_blind_path(setup):
    return setup["paths"]["tier_raw_blind"]


def tier_path(setup):
    return setup["paths"]["tier"]


def get_tier_path(setup, tier):
    if tier in [
        "raw",
        "tcm",
        "dsp",
        "hit",
        "ann",
        "evt",
        "psp",
        "pht",
        "pan",
        "pet",
        "skm",
    ]:
        return setup["paths"][f"tier_{tier}"]
    msg = f"no tier matching:{tier}"
    raise ValueError(msg)


def pars_path(setup):
    return setup["paths"]["par"]


def get_pars_path(setup, tier):
    if tier in ["raw", "tcm", "dsp", "hit", "evt", "psp", "pht", "pet"]:
        return setup["paths"][f"par_{tier}"]
    msg = f"no tier matching:{tier}"
    raise ValueError(msg)


def tmp_par_path(setup):
    return setup["paths"]["tmp_par"]


def tmp_plts_path(setup):
    return setup["paths"]["tmp_plt"]


def plts_path(setup):
    return setup["paths"]["plt"]


def par_overwrite_path(setup):
    return setup["paths"]["par_overwrite"]


def config_path(setup):
    return setup["paths"]["config"]


def chan_map_path(setup):
    return setup["paths"]["chan_map"]


def det_status_path(setup):
    return setup["paths"]["detector_status"]


def metadata_path(setup):
    return setup["paths"]["metadata"]


def detector_db_path(setup):
    return setup["paths"]["detector_db"]


def log_path(setup):
    return setup["paths"]["log"]


def tmp_log_path(setup):
    return setup["paths"]["tmp_log"]


def tmp_benchmark_path(setup):
    return setup["paths"]["tmp_benchmark"]


def filelist_path(setup):
    return setup["paths"]["tmp_filelists"]


# leaf-level paths eligible for symlinking (parent dirs 'tier' and 'par' are
# intentionally excluded to avoid parent/child symlink conflicts)
_SYMLINK_CANDIDATES = (
    "tier_daq",
    "tier_raw",
    "tier_tcm",
    "tier_dsp",
    "tier_hit",
    "tier_ann",
    "tier_evt",
    "tier_psp",
    "tier_pht",
    "tier_pan",
    "tier_pet",
    "tier_skm",
    "tier_raw_blind",
    "par_raw",
    "par_tcm",
    "par_dsp",
    "par_hit",
    "par_evt",
    "par_psp",
    "par_pht",
    "par_pet",
    "plt",
    "metadata",
)


def link_external_paths(
    config,
    workflow_basedir: str | Path,
    *,
    logger: logging.Logger | None = None,
) -> None:
    """Symlink overridden paths back into their canonical default locations.

    When a ``paths.<key>`` entry in ``dataflow-config.yaml`` has been set to a
    location outside the current production tree (e.g. reusing the ``tier_dsp``
    output from another production), this function creates a relative symlink at
    the canonical default location so the ``generated/`` tree maintains a
    consistent layout.

    The canonical defaults are read from the ``dataflow-config.yaml`` found at
    ``<workflow_basedir>/../dataflow-config.yaml``, with ``$_`` substituted by
    the current working directory.  The call is a safe no-op when that file does
    not exist.

    For each candidate key:

    - if the configured path already equals the default, any stale symlink at
      the default location is removed;
    - otherwise a relative symlink is created (or refreshed) at the default
      location pointing to the configured path.

    Real directories at the default location are never modified.
    """
    log_ = logger if logger is not None else logging.getLogger(__name__)
    template = Path(workflow_basedir).parent / "dataflow-config.yaml"
    if not template.is_file():
        return

    default_cfg = yaml.safe_load(template.read_text())
    if not isinstance(default_cfg, dict):
        return

    subst_vars(
        default_cfg,
        var_values={"_": str(Path.cwd().resolve())},
        ignore_missing=True,
    )
    default_paths = default_cfg.get("paths", {})
    config_paths = config.get("paths", {})

    for key in _SYMLINK_CANDIDATES:
        current_str = config_paths.get(key)
        default_str = default_paths.get(key)
        if current_str is None or default_str is None:
            continue

        current = Path(current_str)
        default = Path(default_str)

        # os.path.abspath (not Path.resolve) on purpose: don't follow symlinks
        if os.path.abspath(current) == os.path.abspath(default):  # noqa: PTH100
            if default.is_symlink():
                default.unlink()
            continue

        if default.is_symlink():
            if (default.parent / default.readlink()).resolve() == current.resolve():
                continue
            default.unlink()
        elif default.exists():
            log_.warning(
                "%s is an existing non-symlink path; ignoring override -> %s",
                default,
                current,
            )
            continue

        if not current.exists():
            log_.warning(
                "override target %s does not exist; %s will be a broken symlink",
                current,
                default,
            )
        elif not current.is_dir():
            log_.warning(
                "override target %s is not a directory; linking anyway", current
            )

        rel = Path(os.path.relpath(current, start=default.parent))
        default.parent.mkdir(parents=True, exist_ok=True)
        default.symlink_to(rel, target_is_directory=True)
        log_.info("symlinked %s -> %s", default, rel)
