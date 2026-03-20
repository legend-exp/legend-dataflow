Developer Guide
===============

Repository Structure
--------------------

.. code-block:: text

   legend-dataflow/
   ├── workflow/
   │   ├── Snakefile                  # Main workflow entry point
   │   ├── Snakefile-build-raw        # Separate workflow for raw data building
   │   ├── rules/                     # Snakemake rule files (one per processing stage)
   │   ├── profiles/                  # Execution profiles (default, lngs, sator, ...)
   │   └── src/legenddataflow/
   │       ├── methods/               # Core library: file keys, patterns, calibration grouping
   │       └── scripts/               # Executable scripts called by rules
   │           ├── flow/              # Workflow management scripts
   │           ├── tier/              # Data tier building scripts
   │           └── par/               # Parameter generation scripts
   │               ├── geds/          # HPGe detector parameters
   │               └── spms/          # SiPM detector parameters
   ├── docs/                          # Sphinx documentation
   ├── tests/                         # Test suite
   ├── dataflow-config.yaml           # Default configuration
   └── pyproject.toml                 # Package metadata and dependencies

Snakemake Rules
---------------

The workflow is built around Snakemake rules that define how output files are derived
from input files. Rules are organised into one file per processing stage under
``workflow/rules/``.

Each rule specifies:

- **input** – files or Python callables that resolve input paths from wildcards
- **output** – the file(s) the rule produces
- **params** – additional parameters passed to the script
- **log** – where to write execution logs
- **threads** – number of CPU threads to request
- **script** – the Python script (or shell command) to execute

Full details are in the
`Snakemake rules documentation <https://snakemake.readthedocs.io/en/stable/snakefiles/rules.html>`_.

Rule naming conventions
~~~~~~~~~~~~~~~~~~~~~~~

- ``build_tier_{tier}`` – rules that build a data tier from input data
- ``build_{tier}_pars_{detector}`` – rules that derive calibration parameters for
  a detector type (``geds`` for HPGe, ``spms`` for SiPM)
- Rules operating on partition-level data use the ``psp`` / ``pht`` / ``pan`` / ``pet``
  tier names

Parameter generation pattern
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Most tiers follow the same two-step pattern:

1. **Parameter generation rules** – run on calibration data (``cal`` datatype) to
   derive per-channel calibration parameters. These produce ``par_{tier}.yaml`` files.
2. **Tier building rules** – apply those parameters to both calibration and physics data
   to produce tier files.

For most tiers there are two versions:

- **Run-level** – uses a single calibration run
- **Partition-level** – groups multiple runs together (defined in
  ``cal_groupings.yaml``) for improved statistical precision

Scripts
-------

Scripts live under ``workflow/src/legenddataflow/scripts/`` and are called by rules
via Snakemake's ``script:`` directive. They receive inputs, outputs, and parameters
from Snakemake through the ``snakemake`` object.

Script categories:

- ``flow/`` – workflow management: file discovery, channel lists, run finalisation,
  file merging
- ``tier/`` – data tier building: raw conversion, TCM, event reconstruction, skim
- ``par/geds/`` – HPGe parameter generation scripts organised by tier
  (``raw/``, ``tcm/``, ``psp/``, ``pht/``)
- ``par/spms/`` – SiPM parameter generation scripts

Methods Library
---------------

The ``workflow/src/legenddataflow/methods/`` package provides shared utilities:

- **FileKey** (``FileKey.py``) – parses and generates file keys and patterns from
  wildcard components (experiment, period, run, datatype, timestamp)
- **patterns** (``patterns.py``) – defines file path patterns for each tier and
  processing stage
- **paths** (``paths.py``) – resolves configured paths from the config file
- **CalGrouping** (``cal_grouping.py``) – loads and queries partition-level
  calibration groupings from ``cal_groupings.yaml``
- **ParsKeyResolve** (``create_pars_keylist.py``) – resolves which parameter files
  apply to a given file key, using the parameter validity catalog
- **ParsCatalog** (``pars_loading.py``) – loads and manages parameter catalog files

Adding a new processing stage
------------------------------

To add a new tier ``foo``:

1. **Add path configuration** – add ``tier_foo`` and ``par_foo`` entries to
   ``dataflow-config.yaml`` and the ``paths.py`` helper.

2. **Add file patterns** – add pattern definitions for the new tier in
   ``methods/patterns.py``.

3. **Write the rule file** – create ``workflow/rules/foo.smk`` with:

   - A parameter generation rule (if applicable) that reads calibration data and
     writes ``par_foo.yaml``
   - A tier building rule that applies parameters and writes tier files

4. **Write the script(s)** – add scripts to ``scripts/tier/foo.py`` and/or
   ``scripts/par/geds/foo/`` implementing the processing logic.

5. **Include the rule file** – add ``include: "rules/foo.smk"`` to the main
   ``workflow/Snakefile`` in the appropriate order (after its dependencies).

6. **Update table_format** – add the HDF5 table naming pattern for the new tier to
   the ``table_format`` section of the config.

Adding a new execution environment
------------------------------------

To add a new host or container environment:

1. Add a new key under ``execenv`` in ``dataflow-config.yaml`` with the container
   command, image path, and required environment variables.

2. Add a matching profile directory under ``workflow/profiles/<hostname>/`` with a
   ``config.yaml`` specifying Snakemake options appropriate for that host.

Testing
-------

The test suite uses `pytest <https://docs.pytest.org>`_ and is located in ``tests/``.
Run tests with:

.. code-block:: bash

   pytest tests/

Code style follows PEP 8. The project uses ``ruff`` for linting and formatting, which
is enforced in CI. Run locally with:

.. code-block:: bash

   ruff check .
   ruff format .
