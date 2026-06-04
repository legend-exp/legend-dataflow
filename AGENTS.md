# LEGEND Dataflow — Agent / Developer Reference

## Project Overview

`legend-dataflow` orchestrates the full LEGEND L200 data processing pipeline from raw
digitiser output to physics-ready event data using
[Snakemake](https://snakemake.readthedocs.io/). Hundreds of detector channels are
calibrated and optimised in parallel before physics data are processed.

Full documentation: <https://legend-dataflow.readthedocs.io>

---

## Repository Structure

```text
legend-dataflow/
├── workflow/
│   ├── Snakefile                  # Main workflow entry point
│   ├── Snakefile-build-raw        # Separate workflow for raw data building
│   ├── rules/                     # Snakemake rule files (one per processing stage)
│   ├── profiles/                  # Execution profiles (default, lngs, sator, nersc)
│   └── src/legenddataflow/
│       ├── methods/               # Core library: file keys, patterns, calibration grouping
│       └── scripts/               # Executable scripts called by rules
│           ├── flow/              # File discovery, channel lists, run finalisation
│           ├── tier/              # Data tier building (raw, tcm, evt, skm)
│           └── par/               # Parameter generation scripts
│               ├── geds/          # HPGe detector parameters
│               └── spms/          # SiPM detector parameters
├── docs/                          # Sphinx documentation (RST source)
├── tests/                         # pytest test suite
├── dataflow-config.yaml           # Default site configuration
└── pyproject.toml                 # Package metadata and dependencies
```

Python source lives in `workflow/src/legenddataflow/`; `pyproject.toml` sets
`tool.setuptools.package-dir` accordingly.

---

## Processing Tiers

The pipeline transforms data through successive tiers. For each tier, calibration
parameters are first derived from dedicated calibration runs (`cal` datatype) and then
applied to physics data (`phy` datatype).

| Tier  | Rule file(s)                              | Description                                         |
| ----- | ----------------------------------------- | --------------------------------------------------- |
| `raw` | `rules/raw.smk`                           | Convert DAQ (ORCA/FCIO) to LH5; apply blinding      |
| `tcm` | `rules/tcm.smk`                           | Time Coincidence Map; pulser identification         |
| `dsp` | `rules/dsp.smk`, `dsp_pars_geds/spms.smk` | Digital signal processing; per-channel optimisation |
| `hit` | `rules/hit.smk`, `hit_pars_geds.smk`      | Energy calibration and PSD                          |
| `psp` | `rules/psp.smk`, `psp_pars_geds.smk`      | Partition-level DSP (averaged over multiple runs)   |
| `pht` | `rules/pht.smk`, `pht_pars_geds*.smk`     | Partition-level HIT; QC, A/E, LQ calibrations       |
| `ann` | `rules/ann.smk`                           | ANN cuts for coaxial HPGe detectors                 |
| `evt` | `rules/evt.smk`                           | Event-level reconstruction; cross-talk correction   |
| `skm` | `rules/skm.smk`                           | Final physics skim                                  |

**Dependency graph (simplified):**

```
DAQ files → RAW → TCM → DSP/PSP → HIT/PHT → ANN/PAN → EVT/PET → SKM
```

---

## Key Conventions

### Rule naming

- `build_tier_{tier}` — builds a data tier from input data
- `build_{tier}_pars_{detector}` — derives calibration parameters for a detector type
  (`geds` = HPGe, `spms` = SiPM)
- Partition-level rules use `psp` / `pht` / `pan` / `pet` tier names

### File key structure

Files are identified by components: `experiment-period-run-datatype-timestamp`
(e.g. `l200-p03-r001-phy-20230401T000000Z`). The `FileKey` class in
`methods/FileKey.py` parses and generates these.

### File path patterns

All tier and parameter path patterns are defined in `methods/patterns.py`. New tiers
must add entries here and in `dataflow-config.yaml`.

### Parameter catalogs

Validity YAML files track which parameter set applies to which time range:
`{par_path}/{tier}/validity.yaml`. `ParsKeyResolve` (`methods/create_pars_keylist.py`)
and `ParsCatalog` (`methods/pars_loading.py`) implement catalog resolution.

### Snakemake target format

```
[all|valid|<runlist-key>]-{experiment}-{period}-{run}-{datatype}-{tier}.gen
```

Wildcards (`*`) and multi-value selectors (`_`-separated) work for most components.

---

## Dependencies

**Runtime** (pinned in `pyproject.toml`):

| Package                   | Role                                     |
| ------------------------- | ---------------------------------------- |
| `pygama`                  | Gamma-ray analysis framework             |
| `dspeed`                  | Digital signal processing                |
| `legend-pydataobj`        | LEGEND data objects (LH5 format)         |
| `legend-lh5io`            | LH5 I/O layer                            |
| `legend-daq2lh5`          | DAQ format → LH5 conversion              |
| `legend-dataflow-scripts` | Auxiliary dataflow scripts               |
| `pylegendmeta`            | LEGEND metadata access                   |
| `dbetto`                  | Database-backed parameter store          |
| `snakemake`               | Workflow management (in `runprod` extra) |

**Development extras** (`[dev]` = `[runprod,test]` + `pre-commit`):

```bash
uv pip install -e ".[dev]"
```

**Docs extras** (not installed by default):

```bash
uv pip install -e ".[docs]"
```

Python ≥ 3.11 required; 3.12 is the recommended version.

---

## Code Style

The project uses **ruff** for both linting and formatting (replacing Black/isort).

```bash
ruff check .          # lint
ruff format .         # format
```

Key ruff configuration (from `pyproject.toml`):

- Source root: `workflow/src`
- Enabled rule sets: `ARG`, `B`, `C4`, `EM`, `EXE`, `G`, `I`, `ICN`, `NPY`, `PD`,
  `PGH`, `PIE`, `PL`, `PT`, `PTH`, `RET`, `RUF`, `SIM`, `T20`, `UP`, `YTT`
- `from __future__ import annotations` is required in all Python files
- `T20` (print statements) is allowed in `tests/` and `noxfile.py`

Snakemake rule files are formatted with **snakefmt** (excluded: `channel_merge.smk`).

Shell scripts are checked with **shellcheck**.

---

## Pre-commit

Install hooks after cloning:

```bash
pre-commit install
```

Active hooks (`.pre-commit-config.yaml`):

| Hook                      | Purpose                                        |
| ------------------------- | ---------------------------------------------- |
| `blacken-docs`            | Format Python code blocks in RST/Markdown docs |
| `check-added-large-files` | Block accidental large file commits            |
| `check-yaml/json/toml`    | Validate config file syntax                    |
| `check-merge-conflict`    | Detect unresolved merge markers                |
| `end-of-file-fixer`       | Ensure files end with a newline                |
| `trailing-whitespace`     | Strip trailing spaces                          |
| `name-tests-test`         | Enforce `test_*.py` naming (pytest-first)      |
| `ruff` + `ruff-format`    | Lint and format Python                         |
| `validate-pyproject`      | Validate `pyproject.toml` schema               |
| `check-github-workflows`  | Validate GitHub Actions YAML                   |
| `mypy`                    | Type checking (manual stage only)              |
| `nbstripout`              | Strip notebook outputs before commit           |
| `codespell`               | Spell checking (with physics word exceptions)  |
| `shellcheck`              | Shell script linting                           |
| `rst-backticks` etc.      | RST syntax checks                              |
| `prettier`                | Format YAML, Markdown, JSON, TOML              |
| `snakefmt`                | Format Snakemake rule files                    |

Auto-update schedule: quarterly. `mypy` and `check-manifest` run only on manual stage.

---

## Testing

Tests use **pytest** and live in `tests/`:

```bash
pytest tests/
```

Test files:

- `test_filekey.py` — `FileKey` parsing and pattern generation
- `test_create_pars.py` — parameter key resolution
- `test_pars_loading.py` — parameter catalog loading
- `test_cal_grouping.py` — calibration grouping / partition logic

pytest configuration (from `pyproject.toml`):

- `--strict-markers`, `--strict-config`, `--showlocals`
- `xfail_strict = true` — unexpected passes are failures
- `filterwarnings = ["error"]` — all warnings are errors
- `log_cli_level = "INFO"`

The `tests/dummy_cycle/` and `tests/runprod/` directories contain fixture data for
integration-style tests.

---

## Documentation

Documentation is written in **reStructuredText** under `docs/source/` and built with
**Sphinx** (Furo theme):

```bash
cd docs
make html
```

Source files:

- `index.rst` — top-level table of contents
- `user_manual.rst` — configuration, profiles, running the dataflow
- `pipeline.rst` — processing pipeline overview and tier descriptions
- `developer.rst` — repository structure, conventions, how to extend the pipeline

Published automatically on ReadTheDocs at <https://legend-dataflow.readthedocs.io>.

---

## Extending the Pipeline

### Adding a new DSP processor

If the processor exists in `dspeed`, add it to the relevant config in
`legend-dataflow-config` under `tier/dsp`. Otherwise open a PR to `dspeed` first, or
point `pyproject.toml` at a local version.

### Adding a new calibration script

1. Add a rule to the relevant `.smk` file (e.g. `rules/hit_pars_geds.smk`) writing
   `par_hit_mystep.yaml`.
2. Add the script under `scripts/par/geds/hit/`.

### Adding a new tier `foo`

1. Add `tier_foo` and `par_foo` path entries to `dataflow-config.yaml` and `methods/paths.py`.
2. Add file path patterns in `methods/patterns.py`.
3. Write `workflow/rules/foo.smk` with parameter generation and tier building rules.
4. Write scripts in `scripts/tier/foo.py` and/or `scripts/par/geds/foo/`.
5. Add `include: "rules/foo.smk"` to `workflow/Snakefile` in dependency order.
6. Add the HDF5 table naming pattern to the `table_format` section of the config.

### Adding a new execution environment

1. Add a key under `execenv` in `dataflow-config.yaml` with container command, image
   path, and required environment variables.
2. Add a matching profile under `workflow/profiles/<hostname>/config.yaml`.
