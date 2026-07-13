from __future__ import annotations

import pytest

from legenddataflow.scripts.flow.utils import replace_path
from legenddataflow.scripts.par.geds.pht.util import (
    save_dict_to_files,
    split_files_by_run,
)


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


def test_split_files_by_run_missing_members(tmp_path):
    good = tmp_path / "l200-p03-r000-cal-20230101T000000Z-tier_dsp.lh5"
    good.touch()
    filelist = tmp_path / "cal.filelist"
    filelist.write_text(f"{good}\n{tmp_path}/missing-file.lh5\n")

    with pytest.raises(FileNotFoundError, match="--input-files: input file"):
        split_files_by_run([str(filelist)])

    filelist.write_text(f"{good}\n")
    final_dict, files = split_files_by_run([str(filelist)])
    assert files == [str(good)]
    assert list(final_dict.values()) == [[str(good)]]


def test_save_dict_to_files_creates_parents(tmp_path):
    out = (
        tmp_path
        / "new"
        / "dir"
        / "l200-p03-r000-cal-20230101T000000Z-ch1000000-par_dsp.yaml"
    )
    save_dict_to_files([str(out)], {"20230101T000000Z": {"a": 1}})
    assert out.is_file()
