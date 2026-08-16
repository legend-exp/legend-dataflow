# ruff: noqa: T201, I002

"""Snakemake ``script:`` for the ``gen_filelist`` checkpoint: write the input
files matched by a ``.gen`` target into a filelist."""

from pathlib import Path

try:
    from snakemake.script import snakemake  # snakemake > 8.16
except ImportError:  # not running inside a snakemake job (e.g. unit tests)
    snakemake = None


def write_filelist(files, output_file, label=""):
    print(f"INFO: found {len(files)} files")
    if len(files) == 0:
        print(
            f"WARNING: no files found for the given pattern: {label}. "
            "make sure patterns follows the format: "
            "all-{experiment}-{period}-{run}-{datatype}-{timestamp}-{tier}.gen"
        )
    with Path(output_file).open("w") as f:
        for fn in files:
            f.write(f"{fn}\n")


if snakemake is not None:
    write_filelist(
        list(snakemake.input), snakemake.output[0], label=snakemake.wildcards.label
    )
