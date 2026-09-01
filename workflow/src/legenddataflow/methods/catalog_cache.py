"""In-process cache for dbetto validity catalogs.

Snakemake evaluates ``input:``/``params:`` functions once per job, and several
of them resolve parameter files through a validity catalog. Re-reading the
YAML from disk on every call dominates DAG-build time on large targets, so
catalogs are cached here per path, invalidated when the file's mtime or size
changes.
"""

from __future__ import annotations

from pathlib import Path

from dbetto.catalog import Catalog

_CACHE: dict[str, tuple[tuple[int, int], Catalog]] = {}


def cached_catalog_read(path) -> Catalog:
    """Return the :class:`dbetto.catalog.Catalog` parsed from ``path``, cached
    per path and invalidated when the file's mtime (nanoseconds) or size
    changes.

    The returned catalog is shared between callers: neither it nor the lists
    returned by its ``valid_for()`` may be mutated.
    """
    p = str(path)
    st = Path(p).stat()
    state = (st.st_mtime_ns, st.st_size)
    hit = _CACHE.get(p)
    if hit is None or hit[0] != state:
        hit = (state, Catalog.read_from(p))
        _CACHE[p] = hit
    return hit[1]
