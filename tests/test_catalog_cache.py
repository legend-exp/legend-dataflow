from __future__ import annotations

import os

from dbetto.catalog import Catalog

from legenddataflow.methods import cached_catalog_read

VALIDITY = """\
- valid_from: 20230101T000000Z
  apply:
    - cal/file1.yaml
"""

VALIDITY_UPDATED = """\
- valid_from: 20230101T000000Z
  apply:
    - cal/file1.yaml
    - cal/file2.yaml
"""


def test_cached_catalog_read_returns_same_object(tmp_path):
    validity = tmp_path / "validity.yaml"
    validity.write_text(VALIDITY)

    first = cached_catalog_read(validity)
    assert isinstance(first, Catalog)
    assert first.valid_for("20230201T000000Z") == ["cal/file1.yaml"]
    assert cached_catalog_read(validity) is first
    # str and Path arguments hit the same cache entry
    assert cached_catalog_read(str(validity)) is first


def test_cached_catalog_read_invalidates_on_mtime_change(tmp_path):
    validity = tmp_path / "validity.yaml"
    validity.write_text(VALIDITY)

    first = cached_catalog_read(validity)
    validity.write_text(VALIDITY_UPDATED)
    st = validity.stat()
    os.utime(validity, (st.st_atime, st.st_mtime + 1))

    second = cached_catalog_read(validity)
    assert second is not first
    assert second.valid_for("20230201T000000Z") == [
        "cal/file1.yaml",
        "cal/file2.yaml",
    ]
