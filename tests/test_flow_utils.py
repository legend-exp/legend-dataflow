from __future__ import annotations

from legenddataflow.scripts.flow.utils import replace_path


def test_replace_path_nested():
    d = {
        "file": "/old/dir/file.lh5",
        "nested": {"list": ["/old/dir/file.lh5", 42]},
        "untouched": 1,
    }
    out = replace_path(d, "/old/dir/file.lh5", "/new/dir/merged.lh5")
    assert out["file"] == "$_/merged.lh5"
    assert out["nested"]["list"] == ["$_/merged.lh5", 42]
    assert out["untouched"] == 1


def test_replace_path_by_name():
    # strings containing only the basename of the old path are also replaced
    out = replace_path({"file": "file.lh5"}, "/old/dir/file.lh5", "/new/merged.lh5")
    assert out == {"file": "merged.lh5"}
