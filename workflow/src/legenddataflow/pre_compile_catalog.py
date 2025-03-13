from datetime import datetime, timezone
from pathlib import Path

from dbetto import TextDB
from dbetto.catalog import Catalog


def pre_compile_catalog(validity_path: str | Path):
    if isinstance(validity_path, str):
        validity_path = Path(validity_path)
    catalog = Catalog.read_from(validity_path / "validity.yaml")
    entries = {}
    textdb = TextDB(validity_path, lazy=False)
    for system in catalog.entries:
        entries[system] = []
        for entry in catalog.entries[system]:
            db = textdb.on(
                datetime.fromtimestamp(entry.valid_from, tz=timezone.utc), system=system
            )
            new_entry = Catalog.Entry(entry.valid_from, db)
            entries[system].append(new_entry)
    return Catalog(entries)
