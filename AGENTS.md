# LEGEND Dataflow

## Project Overview

`legend-dataflow` orchestrates the full LEGEND L200 python-based data
processing pipeline from raw digitiser output to physics-ready event data using
[Snakemake](https://snakemake.readthedocs.io/). Hundreds of detector channels
are calibrated and optimised in parallel before physics data are processed.

## Repository Structure

```text
legend-dataflow/
├── workflow/
│   ├── Snakefile                  # Main workflow entry point
│   ├── Snakefile-build-raw        # Separate workflow for raw data building
│   ├── rules/                     # Snakemake rule files (one per processing stage)
│   ├── profiles/                  # Execution profiles (default, lngs, lngs-build-raw, sator)
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

## Processing Tiers

- The output tier data is formatted in LH5 (LEGEND HDF5). Specification
  available [here](https://legend-exp.github.io/legend-data-format-specs/dev/)
- The pipeline transforms data through successive tiers. For each tier,
  calibration parameters are first derived from dedicated calibration runs
  (`cal` datatype) and then applied to physics data (`phy` datatype).

| Tier  | Rule file(s)                              | Description                                         |
| ----- | ----------------------------------------- | --------------------------------------------------- |
| `raw` | `rules/raw.smk`                           | Convert DAQ (ORCA/FCIO) to LH5; apply blinding      |
| `tcm` | `rules/tcm.smk`                           | Time Coincidence Map; pre-compute event structure   |
| `dsp` | `rules/dsp.smk`, `dsp_pars_geds/spms.smk` | Digital signal processing; per-channel optimisation |
| `hit` | `rules/hit.smk`, `hit_pars_geds.smk`      | Energy calibration and PSD                          |
| `psp` | `rules/psp.smk`, `psp_pars_geds.smk`      | Partition-level DSP (averaged over multiple runs)   |
| `pht` | `rules/pht.smk`, `pht_pars_geds*.smk`     | Partition-level HIT; QC, A/E, LQ calibrations       |
| `ann` | `rules/ann.smk`                           | ANN cuts for coaxial HPGe detectors                 |
| `evt` | `rules/evt.smk`                           | Event-level reconstruction; cross-talk correction   |
| `skm` | `rules/skm.smk`                           | Final physics skim                                  |

- Dependency graph (simplified): `DAQ files → RAW → TCM → DSP/PSP → HIT/PHT →
ANN/PAN → EVT/PET → SKM`
- The data can be processed in two modes: by extracting (calibration,
  optimized) parameters from the `cal` data from the same run or from a merger
  of multiple calibration runs (i.e. "partition").

## Official data production

- The data production is usually run on dedicated HPC centres, not on personal
  computers.
- The LNGS center is accessible via SSH host `legend-login.lngs.infn.it`, the
  data productions are stored below `/data2/public/prodenv/`
- The NERSC center is accessible via SSH host `perlmutter-p1.nersc.gov`. The
  data is at `/global/cfs/cdirs/m2676/data/lngs/l200/public/prodenv/`
- You can explore the data produced in different productions there. Ask the
  user about the SSH username to use and how to access relevant software remotely
- There are two categories of productions: `prod-blind` for the blinded data
  and `prod-orig` for the blind data.
- Other categories: `auto`: semi-online automatic basic data production. `tmp`:
  temporary development productions. `ref`: reference productions.

## Key Conventions

- File key structure. Files are identified by components:
  `experiment-period-run-datatype-timestamp` (e.g.
  `l200-p03-r001-phy-20230401T000000Z`). The `FileKey` class in
  `methods/FileKey.py` parses and generates these.
- File path patterns. All tier and parameter path patterns are defined in
  `methods/patterns.py`. New tiers must add entries here and in
  `dataflow-config.yaml`.
- Parameter catalogs. Validity YAML files track which parameter set applies to
  which time range: `{par_path}/{tier}/validity.yaml`. `ParsKeyResolve`
  (`methods/create_pars_keylist.py`) and `ParsCatalog`
  (`methods/pars_loading.py`) implement catalog resolution.
- Snakemake target format:
  `[all|valid|<runlist-key>]-{experiment}-{period}-{run}-{datatype}-{tier}.gen`
  Wildcards (`*`) and multi-value selectors (`_`-separated) work for most
  components.
- Rule naming. `build_tier_{tier}` — builds a data tier from input data.
  `build_{tier}_pars_{detector}` — derives calibration parameters for a
  detector type (`geds` = HPGe, `spms` = SiPM). Partition-level rules use `psp`
  / `pht` / `pan` / `pet` tier names

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

- `uv pip install .[runprod]` to get everything needed to run the production
- `dataprod -v install dataflow-config.yaml` to install the software in the
  production environment
- `uv pip install .[dev]` to get everything needed to develop

## Metadata

- This is the configuration of the dataflow and other information about the
  data
- Stored in the legend-exp/legend-metadata git repository, usually cloned below
  `./inputs` (path configured in the dataflow config)
- Refer to agent files there and subfolders for more details
- Typically useful to configure an AI agent specialized to the metadata

## Linting

- run `pre-commit install` always before any commit and address any flagged
  issue.

## Testing

- Tests use **pytest** and live in `tests/`
- The `tests/dummy_cycle/` and `tests/runprod/` directories contain fixture
  data for integration-style tests.

## Documentation

- Get dependencies: `uv pip install -e ".[docs]"`
- Documentation is written in **reStructuredText** under `docs/source/` and
  built with **Sphinx**: `cd docs && make`

## Boundaries

- ALWAYS run pre-commit before committing
- Do not make statements about features of the LEGEND data unless you have
  verified them on data
- Make sure changes to tier processing are propagated to both standard and
  partition tiers
