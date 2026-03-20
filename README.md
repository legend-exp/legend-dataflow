# LEGEND L200 dataflow

Implementation of an automatic data processing pipeline for
[LEGEND](https://www.legend-experiment.org/) (Large Enriched Germanium Experiment for
Neutrinoless double beta decay) L200 data, based on
[Snakemake](https://snakemake.readthedocs.io/).

## Overview

`legend-dataflow` orchestrates the full chain of data processing from raw digitiser
output to physics-ready event data. It calibrates and optimises hundreds of detector
channels in parallel before bringing them together to process physics data.

The pipeline processes data through a series of tiers:

| Tier  | Description                                           |
| ----- | ----------------------------------------------------- |
| `raw` | Converted from DAQ format (ORCA/FCIO) to LH5          |
| `tcm` | Timing correction module                              |
| `dsp` | Digital signal processing (pole-zero, energy filters) |
| `hit` | Energy calibration and pulse shape discrimination     |
| `psp` | Partition-level DSP (averaged over calibration runs)  |
| `pht` | Partition-level HIT with advanced quality control     |
| `ann` | Artificial neural network cuts (coax detectors only)  |
| `evt` | Event-level reconstruction and cross-talk correction  |
| `skm` | Final physics skim                                    |

For each tier, calibration parameters are first derived per-channel from dedicated
calibration data, then applied to physics data in parallel.

## Installation

Clone the repository and set up a virtual environment (requires Python 3.11+):

```bash
git clone https://github.com/legend-exp/legend-dataflow.git
cd legend-dataflow
uv venv --python 3.12
source .venv/bin/activate
uv pip install -e ".[dev]"
```

Adapt `dataflow-config.yaml` to your site (paths, execution environment), then install
the software environment:

```bash
dataflow -v install -s <host> dataflow-config.yaml
```

where `<host>` is one of the execution environments defined in the config (`bare`,
`lngs`, `sator`, `nersc`).

## Running the dataflow

The main target format is:

```
[all|valid]-{experiment}-{period}-{run}-{datatype}-{tier}.gen
```

Examples:

```bash
# Process all physics data through the SKM tier
snakemake --profile workflow/profiles/default all-l200-p03-r001-phy-skm.gen

# Process a specific period and datatype through DSP
snakemake --profile workflow/profiles/default sel-l200-p03-*-cal-dsp.gen

# Process multiple runs
snakemake --profile workflow/profiles/default all-l200-p03-r000_r001-phy-hit.gen
```

Use `all` to process all data or `valid` to process only analysis-selected data 
(any keyword present in the `runlists.yaml` file in `legend-datasets` can be used here).
Wildcards (`*`) and multi-value selectors (`_`-separated) are supported for most
label components.

## Documentation

Full documentation is available at the project's ReadTheDocs page:

- **User Manual** – configuration, profiles, and running the dataflow
- **Pipeline Overview** – how the processing pipeline works step by step
- **Developer Guide** – how to extend and contribute to the pipeline
- **API Reference** – Python API documentation

## Related projects

- [legend-pydataobj](https://legend-pydataobj.readthedocs.io) – LEGEND data objects (LH5 format)
- [legend-daq2lh5](https://legend-daq2lh5.readthedocs.io) – DAQ format conversion to LH5
- [dspeed](https://dspeed.readthedocs.io) – Digital signal processing
- [pygama](https://pygama.readthedocs.io) – Gamma-ray analysis

## License

MIT License. See [LICENSE.md](LICENSE.md).
