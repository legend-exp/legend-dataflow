User Manual
===========

Installation
============

Prerequisites
-------------

- Python 3.11 or later
- `uv <https://docs.astral.sh/uv/>`_ (or another virtual environment manager)
- Access to the LEGEND metadata repository

Clone and install the package:

.. code-block:: bash

   git clone https://github.com/legend-exp/legend-dataflow.git
   cd legend-dataflow
   uv venv --python 3.12
   source .venv/bin/activate
   uv pip install -e ".[dev]"

The ``[dev]`` extras include development tools such as testing and linting
dependencies. For production use, omit ``[dev]``.

Installing the software environment
------------------------------------

After configuring ``dataflow-config.yaml`` for your site (see :ref:`configuration`
below), install the execution environment:

.. code-block:: bash

   dataflow -v install -s <host> dataflow-config.yaml

where ``<host>`` is one of the execution environments defined in the config (e.g.
``bare``, ``lngs``, ``sator``, ``nersc``). This installs all required software into
``.snakemake/legend-dataflow/venv``.

.. note::

   If you update the software, clear the numba cache directory (defined in the config
   under ``execenv.<host>.env.NUMBA_CACHE_DIR``) to avoid stale compiled code.


.. _configuration:

Configuration
=============

Data processing resources are configured via a single site-dependent YAML file,
conventionally named ``dataflow-config.yaml``. The default file in the repository root
can serve as a starting point.

Key settings
------------

+-------------------------------+--------------------------------------------------------------+
| Parameter                     | Description                                                  |
+===============================+==============================================================+
| ``legend_metadata_version``   | Version tag of legend-metadata to check out automatically.  |
|                               | Set a custom path at ``paths/metadata`` to use a local copy.|
+-------------------------------+--------------------------------------------------------------+
| ``allow_none_par``            | If ``false``, the workflow aborts when calibration parameter |
|                               | generation fails. If ``true``, it continues with default or  |
|                               | overridden parameters.                                       |
+-------------------------------+--------------------------------------------------------------+
| ``build_file_dbs``            | Whether to generate ``pygama.flow.FileDB`` databases after   |
|                               | each successful production run.                              |
+-------------------------------+--------------------------------------------------------------+
| ``check_log_files``           | Whether to scan log files for errors/warnings after each run.|
+-------------------------------+--------------------------------------------------------------+
| ``multiprocess``              | Enable parallel processing within a single Snakemake job.    |
+-------------------------------+--------------------------------------------------------------+

Paths
-----

All path values support the ``$_`` placeholder, which is substituted with the value of
the ``$PRODENV`` environment variable at runtime. This allows a single config file to
be used across different machines by setting ``PRODENV`` appropriately.

The key path categories are:

- **Input paths**: ``metadata``, ``config``, ``par_overwrite``, ``chan_map``,
  ``detector_status``, ``detector_db``
- **Output tier paths**: ``tier_raw``, ``tier_tcm``, ``tier_dsp``, ``tier_hit``,
  ``tier_psp``, ``tier_pht``, ``tier_ann``, ``tier_evt``, ``tier_pan``, ``tier_pet``,
  ``tier_skm``
- **Parameter directories**: ``par_raw``, ``par_tcm``, ``par_dsp``, ``par_hit``,
  ``par_psp``, ``par_pht``, ``par_pet``
- **Scratch/log directories**: ``tmp_plt``, ``tmp_log``, ``tmp_filelists``, ``tmp_par``,
  ``log``, ``plt``

Execution environments
----------------------

The ``execenv`` section defines how scripts are executed. Each named environment
specifies an optional container command and the environment variables to set:

- ``bare`` – run directly in the local Python environment (no container)
- ``lngs`` – Apptainer container on the LNGS cluster
- ``sator`` – Apptainer container on the Sator cluster
- ``nersc`` – Shifter container at NERSC

To add a new environment, add a new key under ``execenv`` following the same pattern.


Profiles
========

Snakemake execution profiles are stored in ``workflow/profiles/``. Each profile is a
directory containing a ``config.yaml`` that sets Snakemake options such as the number
of cores and memory constraints.

The available profiles are:

- ``default`` – bare-metal execution using all available cores
- ``lngs`` – LNGS computing cluster settings
- ``sator`` – Sator computing cluster settings
- ``lngs-build-raw`` – settings specific to raw data building at LNGS

Specify a profile with the ``--profile`` flag:

.. code-block:: bash

   snakemake --profile workflow/profiles/lngs all-l200-p03-r001-phy-skm.gen

A full list of configurable Snakemake options is available in the
`Snakemake CLI documentation <https://snakemake.readthedocs.io/en/stable/executing/cli.html>`_.


Running the Dataflow
====================

The ``$PRODENV`` environment variable must be set to the root of your production
environment before running Snakemake:

.. code-block:: bash

   export PRODENV=/path/to/your/production/environment

Single-file targets
-------------------

At the most basic level you can ask Snakemake to build a single output file, and it
will work out all the dependencies automatically:

.. code-block:: bash

   snakemake /path/to/generated/tier/dsp/p03/r000/l200-p03-r000-cal-20230401T000000Z-tier_dsp.lh5

Batch targets with ``.gen`` files
-----------------------------------

In practice, you will want to process many files at once. The special ``.gen`` target
format triggers processing of all matching files up to a given tier.

The target format is:

.. code-block:: text

   [all|valid]-{experiment}-{period}-{run}-{datatype}-{tier}.gen

where:

- ``all`` / ``valid`` – process all data, or only data selected for analysis, 
any keyword in `runlists.yaml` in `legend-datasets` is a possible option.
- ``experiment`` – experiment name (e.g. ``l200``)
- ``period`` – data-taking period (e.g. ``p03``)
- ``run`` – run number (e.g. ``r001``)
- ``datatype`` – data type (e.g. ``phy`` for physics, ``cal`` for calibration)
- ``tier`` – output tier to build up to (e.g. ``dsp``, ``hit``, ``skm``)

Any component except ``tier`` can be replaced by a wildcard (``*``) to match all
values, or a ``_``-separated list to match multiple specific values.

Examples:

.. code-block:: bash

   # Process all physics data from period p03, run r001, through to the SKM tier
   snakemake all-l200-p03-r001-phy-skm.gen

   # Process all calibration data from any run in period p03 to the DSP tier
   snakemake all-l200-p03-*-cal-dsp.gen

   # Process physics data from runs r000 and r001 to the HIT tier
   snakemake all-l200-p03-r000_r001-phy-hit.gen

   # Process analysis-selected physics data from any period and run to SKM
   snakemake valid-l200-*-*-phy-skm.gen

On success, the empty marker file ``{label}-{tier}.gen`` is created to record
that production completed successfully.

Post-processing
---------------

On successful completion, the workflow automatically:

- Collects warnings and errors from individual log files into a summary log
- Generates a Snakemake HTML report saved under the log directory
- Builds ``pygama.flow.FileDB`` databases for the output files (if enabled)
- Generates lists of valid file keys
- Writes parameter validity catalog files for each tier

Monitoring
==========

You can use the snkmt TUI for monitoring. Available with `snkmt --console`


Software Containers
===================

The dataflow uses container environments for reproducible execution on HPC systems.
Rather than Snakemake's built-in Singularity/Apptainer support, it manages containers
through the ``execenv`` configuration, giving finer control over which commands are
containerised.

Container settings are only required when Snakemake itself runs outside the container
(e.g. when submitting jobs to a batch system). If the entire workflow runs inside the
container, no special container configuration is needed.

Supported container runtimes:

- **Apptainer** (formerly Singularity) – used at LNGS and Sator
- **Shifter** – used at NERSC

The container image path and runtime command are configured per environment in the
``execenv`` section of ``dataflow-config.yaml``.


Parameter Overrides
===================

Calibration parameters can be overridden on a per-channel, per-run basis by placing
override files in the directory specified by ``paths/par_overwrite`` in the
configuration. These take precedence over parameters derived by the automatic
calibration pipeline.

The override directory follows the same hierarchical structure as the parameter output
directories.


Run Validity and Ignored Cycles
=================================

Two files in the ``detector_status`` directory control which data are included:

- ``ignored_daq_cycles.yaml`` – lists DAQ cycles to exclude from processing entirely
- ``run_override.yaml`` – allows overriding the validity window for specific runs,
e.g. to apply a previous valid calibration to a subsequent run

These files are part of the ``legend-datasets`` repository.
