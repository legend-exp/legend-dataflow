from __future__ import annotations

import subprocess
from pathlib import Path

from legenddataflow.methods import metadata_updates as mu
from legenddataflow.methods.metadata_updates import check_metadata_updates


def git(repo, *args, check=True):
    r = subprocess.run(
        ["git", "-C", str(repo), *args], capture_output=True, text=True, check=False
    )
    if check and r.returncode:
        msg = f"git {' '.join(args)} failed in {repo}: {r.stderr}"
        raise RuntimeError(msg)
    return r


def init_repo(path: Path):
    path.mkdir(parents=True, exist_ok=True)
    git(path, "init", "-q", "-b", "main")
    git(path, "config", "user.email", "test@example.com")
    git(path, "config", "user.name", "test")
    git(path, "config", "commit.gpgsign", "false")


def commit_file(repo, rel, content, msg):
    p = Path(repo) / rel
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(content)
    git(repo, "add", "-A")
    git(repo, "commit", "-q", "-m", msg)


class CapLogger:
    """Minimal logger capturing .info/.warning lines (with %-formatting)."""

    def __init__(self):
        self.lines = []

    def info(self, m, *a):
        self.lines.append("INFO: " + (m % a if a else m))

    def warning(self, m, *a):
        self.lines.append("WARN: " + (m % a if a else m))


def build(base: Path, branch: str = "main"):
    """Create an ``upstream`` repo (main + optional branch) and a ``work`` clone on ``branch``."""
    up = base / "upstream"
    init_repo(up)
    commit_file(up, "statuses/seed.yaml", "seed: 1\n", "seed")
    commit_file(up, "raw/validity.yaml", "- valid_from: 0\n", "seed raw")
    if branch != "main":
        git(up, "branch", branch)
    work = base / "work"
    git(base, "clone", "-q", str(up), str(work))
    git(work, "config", "user.email", "test@example.com")
    git(work, "config", "user.name", "test")
    git(work, "config", "commit.gpgsign", "false")
    git(work, "checkout", "-q", branch)
    return up, work


def cfg(repo):
    # Point both metadata targets at the one scratch repo; we assert on the
    # "legend-datasets" label (its subtree is statuses/, which our commits touch).
    return {"paths": {"detector_status": str(repo), "par_overwrite": str(repo)}}


def test_clean_behind_fast_forwards(tmp_path):
    up, work = build(tmp_path, "main")
    commit_file(
        up, "statuses/new.yaml", "x: 1\n", "upstream advance"
    )  # main moves ahead
    log = CapLogger()
    res = check_metadata_updates(cfg(work), mode="ff-only", logger=log)[
        "legend-datasets"
    ]
    assert res["updated"] is True
    assert res["behind_branch"] == 0
    assert res["dirty"] is False


def test_dirty_behind_warns_only(tmp_path):
    up, work = build(tmp_path, "main")
    commit_file(up, "statuses/new.yaml", "x: 1\n", "upstream advance")
    (work / "statuses/seed.yaml").write_text("seed: LOCAL EDIT\n")  # make work dirty
    log = CapLogger()
    res = check_metadata_updates(cfg(work), mode="ff-only", logger=log)[
        "legend-datasets"
    ]
    assert res["updated"] is False
    assert res["dirty"] is True
    assert res["behind_branch"] >= 1
    assert any("DIRTY" in line for line in log.lines)


def test_fix_on_main_flagged_not_merged(tmp_path):
    up, work = build(tmp_path, "p16")  # work tracks period branch p16
    # a fix lands on upstream main, touching statuses/ -- not on p16
    git(up, "checkout", "-q", "main")
    commit_file(
        up, "statuses/fix.yaml", "V06649M: {processable: false}\n", "fix on main"
    )
    log = CapLogger()
    res = check_metadata_updates(cfg(work), mode="ff-only", logger=log)[
        "legend-datasets"
    ]
    assert res["branch"] == "p16"
    assert res["behind_branch"] == 0
    assert res["behind_main"] >= 1
    assert res["updated"] is False


def test_detached_head_skipped(tmp_path):
    _up, work = build(tmp_path, "main")
    head = git(work, "rev-parse", "HEAD").stdout.strip()
    git(work, "checkout", "-q", head)  # detached HEAD
    log = CapLogger()
    res = check_metadata_updates(cfg(work), mode="ff-only", logger=log)[
        "legend-datasets"
    ]
    assert res["updated"] is False
    assert res["branch"] in ("", "HEAD")
    assert any("detached" in line for line in log.lines)


def test_offline_fetch_failure_skipped(tmp_path):
    _up, work = build(tmp_path, "main")
    git(
        work, "remote", "set-url", "origin", str(tmp_path / "does-not-exist")
    )  # break fetch
    log = CapLogger()
    res = check_metadata_updates(cfg(work), mode="ff-only", logger=log)[
        "legend-datasets"
    ]
    assert res["updated"] is False
    assert any("fetch failed" in line for line in log.lines)


def test_missing_repo_skipped(tmp_path):
    log = CapLogger()
    res = check_metadata_updates(cfg(tmp_path / "nope"), mode="ff-only", logger=log)[
        "legend-datasets"
    ]
    assert res["updated"] is False
    assert any("no git repo" in line for line in log.lines)


def test_git_helper_never_raises_on_missing_binary(tmp_path, monkeypatch):
    # A missing git executable must come back as a non-zero CompletedProcess,
    # not an OSError, so the onstart hook can warn-and-skip.
    real_run = subprocess.run

    def fake_run(cmd, *a, **k):
        cmd = ["/nonexistent/git-binary", *cmd[1:]]
        return real_run(cmd, *a, **k)

    monkeypatch.setattr(mu.subprocess, "run", fake_run)
    cp = mu._git(tmp_path, "status")
    assert cp.returncode != 0
