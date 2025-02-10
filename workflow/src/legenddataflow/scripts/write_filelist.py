# ruff: noqa: F821, T201
from pathlib import Path

from snakemake.script import snakemake  # snakemake > 8.16

print(f"INFO: found {len(snakemake.input)} files")
if len(snakemake.input) == 0:
    print(
        f"WARNING: no files found for the given pattern: {snakemake.wildcards.label}. "
        "make sure patterns follows the format: "
        "all-{experiment}-{period}-{run}-{datatype}-{timestamp}-{tier}.gen"
    )
with Path(snakemake.output[0]).open("w") as f:
    for fn in snakemake.input:
        f.write(f"{fn}\n")
