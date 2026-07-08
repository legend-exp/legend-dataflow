# Contributing

Contributions are welcome via pull requests against `main`.

## Development setup

```bash
git clone https://github.com/legend-exp/legend-dataflow.git
cd legend-dataflow
uv venv --python 3.12
source .venv/bin/activate
uv pip install -e ".[dev]"
pre-commit install
```

## Code style

Linting and formatting are enforced by [pre-commit](https://pre-commit.com/):
ruff (lint + format) for Python, snakefmt for Snakemake files, shellcheck for
shell scripts, and prettier for YAML/Markdown/JSON/TOML. pre-commit.ci runs on
every pull request and pushes autofix commits, but running the hooks locally
gives faster feedback:

```bash
pre-commit run --all-files
```

Docstrings use the NumPy style — the API reference is auto-generated from them
with Sphinx, so undocumented code shows up as empty stubs in the published docs.

Commit messages loosely follow conventional-commit prefixes, e.g.
`fix(evt): ...`, `docs: ...`, `chore: ...`.

## Tests

```bash
python -m pytest              # unit tests
./tests/runprod/run-all.sh    # integration tests (from the repo root;
                              # requires LEGEND_METADATA to be set)
```

The integration tests validate the Snakemake DAG wiring with dummy files; real
data processing (~1 TB per run) is far beyond CI capacity, so new processing
logic should be covered by unit tests where practical.

## Documentation

```bash
uv pip install -e ".[docs]"
make -C docs                  # builds HTML into docs/build (warnings are errors)
```

See the [Developer Guide](https://legend-dataflow.readthedocs.io) for
step-by-step instructions on extending the pipeline (new processors,
calibration scripts, tiers, and execution environments).

## AI-assisted contributions

Permitted with human review and disclosure — see [AI_POLICY.md](AI_POLICY.md).
