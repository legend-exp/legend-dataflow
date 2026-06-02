"""
This module contains resolvers for the config.json dictionary
"""

from __future__ import annotations

import logging
import os
from pathlib import Path


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


# Maps each paths.<key> candidate to its canonical default location relative to
# the production root (CWD).  Parent keys 'tier' and 'par' are intentionally
# excluded to avoid parent/child symlink conflicts.
_DEFAULT_SUBPATHS: dict[str, str] = {
    "tier_daq": "generated/tier/daq",
    "tier_raw": "generated/tier/raw",
    "tier_tcm": "generated/tier/tcm",
    "tier_dsp": "generated/tier/dsp",
    "tier_hit": "generated/tier/hit",
    "tier_ann": "generated/tier/ann",
    "tier_evt": "generated/tier/evt",
    "tier_psp": "generated/tier/psp",
    "tier_pht": "generated/tier/pht",
    "tier_pan": "generated/tier/pan",
    "tier_pet": "generated/tier/pet",
    "tier_skm": "generated/tier/skm",
    "tier_raw_blind": "generated/tier/raw-blind",
    "par_raw": "generated/par/raw",
    "par_tcm": "generated/par/tcm",
    "par_dsp": "generated/par/dsp",
    "par_hit": "generated/par/hit",
    "par_evt": "generated/par/evt",
    "par_psp": "generated/par/psp",
    "par_pht": "generated/par/pht",
    "par_pet": "generated/par/pet",
    "plt": "generated/plt",
    "metadata": "inputs",
}


def link_external_paths(
    config,
    *,
    logger: logging.Logger | None = None,
) -> None:
    """Symlink overridden paths back into their canonical default locations.

    When a ``paths.<key>`` entry in ``dataflow-config.yaml`` has been set to a
    location outside the current production tree (e.g. reusing the ``tier_dsp``
    output from another production), this function creates a relative symlink at
    the canonical default location so the ``generated/`` tree maintains a
    consistent layout.

    The canonical default for each key is derived from :data:`_DEFAULT_SUBPATHS`
    relative to the current working directory (the production root).

    For each candidate key:

    - if the configured path already equals the default, any stale symlink at
      the default location is removed;
    - otherwise a relative symlink is created (or refreshed) at the default
      location pointing to the configured path.

    Real directories at the default location are never modified.
    """
    log_ = logger if logger is not None else logging.getLogger(__name__)
    config_paths = config.get("paths", {})
    root = Path.cwd()

    for key, subpath in _DEFAULT_SUBPATHS.items():
        current_str = config_paths.get(key)
        if current_str is None:
            continue

        current = Path(current_str)
        default = root / subpath

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
        log_.info("symlinking %s -> %s", default, rel)
        default.parent.mkdir(parents=True, exist_ok=True)
        default.symlink_to(rel, target_is_directory=True)
