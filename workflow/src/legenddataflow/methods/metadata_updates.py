"""Check metadata submodules for upstream updates at workflow start.

Invoked from the ``onstart`` hooks of ``workflow/Snakefile`` and
``workflow/Snakefile-build-raw`` (guarded by the ``auto_update_metadata`` config
knob and skipped under ``--touch``). See :func:`check_metadata_updates`.

Production advances each metadata submodule on its own branch, in place
(legend-datasets -> a period branch e.g. p16/p17; legend-dataflow-overrides ->
main/new_auto), so a merged fix only reaches a run when that submodule is pulled.
At workflow start we fetch (read-only) and, per submodule:
  * fast-forward it iff it is CLEAN and strictly behind its tracked branch;
  * otherwise only warn -- never clobber the dirty/ahead/detached working trees
    production keeps (e.g. datasets' local un-pushed edits);
  * also flag commits on origin/main (touching the relevant subtree) that are not
    yet on the tracked branch -- e.g. a legend-datasets fix merged to main but not
    forward-ported onto the period branch.

The only mutation is the opt-in fast-forward merge. ``fetch`` is time-bounded and
failure-tolerant (a timeout, a missing ``git`` binary, or an offline / no-permission
run simply skips the check), so the hook never aborts the workflow.

NOTE: metadata is loaded at parse time (top of the Snakefile), before ``onstart``
runs, so a fast-forward here takes effect on the NEXT run, not the current one.
For the hourly production this is fine (a fix lands within an hour) and avoids
changing metadata mid-parse.
"""

from __future__ import annotations

import logging
import subprocess
from pathlib import Path

log = logging.getLogger(__name__)


def _git(repo: Path, *args: str, timeout: int = 120) -> subprocess.CompletedProcess:
    """Run ``git -C <repo> <args>``, capturing output; never raises.

    A non-zero exit, a timeout, or a missing ``git`` executable all come back as
    a :class:`subprocess.CompletedProcess` with a non-zero ``returncode``, so
    callers can uniformly warn-and-skip instead of letting the onstart hook abort
    the workflow.
    """
    cmd = ["git", "-C", str(repo), *args]
    try:
        return subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=timeout,
            check=False,
        )
    except (subprocess.TimeoutExpired, OSError) as exc:
        return subprocess.CompletedProcess(
            cmd, returncode=1, stdout="", stderr=str(exc)
        )


def _behind(repo: Path, ref: str, subtree: tuple[str, ...] = ()) -> int:
    """Number of commits in ``HEAD..<ref>`` (optionally restricted to a subtree)."""
    args = ["rev-list", "--count", f"HEAD..{ref}"]
    if subtree:
        args += ["--", *subtree]
    out = _git(repo, *args).stdout.strip()
    return int(out) if out.isdigit() else 0


def check_metadata_updates(
    config, *, mode: str = "ff-only", logger=None, notify=None
) -> dict:
    """Fetch + report/refresh the metadata submodules. Returns ``{label: info}``.

    Parameters
    ----------
    config
        The workflow config (needs ``config["paths"]["detector_status"]`` and
        ``config["paths"]["par_overwrite"]``).
    mode
        ``"ff-only"`` fast-forwards clean submodules that are strictly behind;
        ``"warn"`` reports only, never modifies.
    logger
        Object with ``.info`` / ``.warning`` (Snakemake's global ``logger``); if
        ``None``, falls back to this module's :mod:`logging` logger.
    notify
        Optional callable ``notify(str)``, e.g. to post to a Slack webhook.
    """
    out = logger if logger is not None else log
    info_ = out.info
    warn_ = out.warning

    # (label, repo path, subtree that matters for the origin/main divergence check)
    targets = [
        (
            "legend-datasets",
            config["paths"]["detector_status"],
            ("statuses", "run_override.yaml", "ignored_daq_cycles.yaml"),
        ),
        (
            "legend-dataflow-overrides/raw",
            config["paths"]["par_overwrite"],
            ("raw",),
        ),
    ]

    results: dict = {}
    for label, repo_path, subtree in targets:
        repo = Path(repo_path)
        info = {
            "branch": None,
            "behind_branch": 0,
            "behind_main": 0,
            "dirty": False,
            "updated": False,
        }

        if not (repo / ".git").exists():
            warn_("metadata %s: no git repo at %s -- skipping", label, repo)
            results[label] = info
            continue

        branch = _git(repo, "rev-parse", "--abbrev-ref", "HEAD").stdout.strip()
        info["branch"] = branch
        info["dirty"] = bool(_git(repo, "status", "--porcelain").stdout.strip())

        # Read-only fetch of the tracked branch (+ main for the divergence check).
        refs = ([branch] if branch not in ("", "HEAD") else []) + ["main"]
        if _git(repo, "fetch", "--quiet", "origin", *refs).returncode != 0:
            warn_(
                "metadata %s: fetch failed (offline / no permission) -- skipping", label
            )
            results[label] = info
            continue

        if branch in ("", "HEAD"):
            warn_(
                "metadata %s: detached HEAD -- cannot track a branch; skipping", label
            )
            results[label] = info
            continue

        info["behind_branch"] = _behind(repo, f"origin/{branch}")
        info["behind_main"] = (
            0 if branch == "main" else _behind(repo, "origin/main", subtree)
        )

        # Fast-forward only if the tree is clean and strictly behind its own branch.
        if mode == "ff-only" and info["behind_branch"] and not info["dirty"]:
            if _git(repo, "merge", "--ff-only", f"origin/{branch}").returncode == 0:
                _git(repo, "submodule", "update", "--init", "--recursive")
                info["updated"] = True
                info["behind_branch"] = 0
            else:
                warn_(
                    "metadata %s: cannot fast-forward %s (diverged) -- left as-is",
                    label,
                    branch,
                )

        if info["behind_branch"] or info["behind_main"] or info["dirty"]:
            msg = (
                f"metadata {label} on {branch}: {info['behind_branch']} behind origin/{branch}"
                + (" (fast-forwarded)" if info["updated"] else "")
                + (
                    f"; {info['behind_main']} fix commit(s) on origin/main not yet on {branch}"
                    if info["behind_main"]
                    else ""
                )
                + (" [DIRTY local edits -- not touched]" if info["dirty"] else "")
            )
            warn_(msg)
            if notify:
                notify(msg)
        else:
            info_("metadata %s on %s: up to date", label, branch)

        results[label] = info

    return results
