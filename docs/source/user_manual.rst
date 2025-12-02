User Manual
-----------

Configuration
=============

Data processing resources are configured via a single site-dependent (and
possibly user-dependent) configuration file, generally named ``dataflow-config.json``.
Although you can choose any arbitrary name.
Edit this file and adjust paths adjusted as necessary. Note that, when running Snakemake,
the default path to the config file is ``./dataflow-config.json``.

The following (non-exhaustive) table shows a list of options:

+---------------------------+--------------------------------------------------------------------------+
| Parameter                 | Description                                                              |
+===========================+==========================================================================+
| legend_metadata_version   | The version of legend_metadata to be used automatically                  |
|                           | (use custom legend_metadata: put one at the ``paths/metadata`` location) |
+---------------------------+--------------------------------------------------------------------------+
| allow_none_par            | if pargen should be run                                                  |
+---------------------------+--------------------------------------------------------------------------+
| paths                     | Paths to legend_metadata, data input and output.                         |
|                           | Adapt e.g. ``paths/raw`` to point to existing raw input data.            |
+---------------------------+--------------------------------------------------------------------------+

Profiles
========

A number of profiles are also included in the ``workflow/profiles`` directory. If none
are specified, the default profile is used. The profile can be specified by
using the ``--profile`` option when running Snakemake. These control how many
jobs are run simultaneously, based on how many cores are specified and the
memory constraints of the system. A full list of all the options that can be
specified to snakemake can be found at `snakemake
<https://snakemake.readthedocs.io/en/stable/executing/cli.html>`_.


Running the Dataflow
====================

To run the dataflow at the most basic level all that is necassary is to tell
snakemake the target file generation. In a simple case this may just be a
single file e.g.
```shell
$ snakemake /data2/public/prodenv/prod-blind/ref-v1.0.0/generated/tier/dsp/p03/r000/l200-p03-r000-cal-20230401T000000Z-tier_dsp.lh5
```
This would generate the file and all the files that are required to generate
it.  In general though we want to generate a large number of files, and we can
do this using the ``gen`` target.

Main output generation
======================

Usually, the main output will be determined by a file-list.  The special output
target ``{label}-{tier}.gen`` is used to generate all files that follow the
label up to the specified tier.  The label is composed of the following parts:

- the filelist designator: in most cases this will be ``all``, but other
  options are specified in the ``runlists.yaml`` file in the `legend-datasets
  <https://github.com/legend-exp/legend-datasets>`_ repository.
- experiment: the experiment name i.e. l200
- period: the period of the data e.g. p03
- run: the run number e.g. r000
- datatype: the data type e.g. cal
- timestamp: the timestamp of the data e.g. 20230401T000000Z

Example:

```shell
$ snakemake all-l200-p03-r001-cal-20230401T000000Z-dsp.gen
```

You can specify as many or as few of these as they like e.g.
``all-l200-p03-dsp.gen`` If you want to specify a lower part of the label but
leave a higher part free, you can use the ``*``` character e.g.
``all-l200-p03-*-cal-dsp.gen`` .  Additionally if you want to specify multiple
options for a part of the label you can use the ``_`` character between e.g.
``all-l200-p03-r000_r001-dsp.gen``.

After the files are created, the empty file ``{label}-{tier}.gen```` will be
created to mark the successful data production.


Monitoring
==========

Snakemake supports monitoring by connecting to a
`panoptes <https://github.com/panoptes-organization/panoptes>`_ server.

Run (e.g.)
```shell
$ panoptes --port 5000
```
in the background to run a panoptes server instance, which comes with a
GUI that can be accessed with a web-brower on the specified port.

Then use the Snakemake option ``--wms-monitor`` to instruct Snakemake to push
progress information to the panoptes server:
```shell
snakemake --wms-monitor http://127.0.0.1:5000 [...]
```

Using software containers
=========================

This dataflow doesn't use Snakemake's internal Singularity support, but
instead supports Singularity containers via
`venv <https://github.com/oschulz/singularity-venv>`_ environments
for greater control.

To use this, the path to ``venv`` and the name of the environment must be set
in ``config.json``.

This is only relevant then running Snakemake *outside* of the software
container, e.g. when using a batch system (see below). If Snakemake
and the whole workflow is run inside of a container instance, no
container-related settings in ``config.json`` are required.
