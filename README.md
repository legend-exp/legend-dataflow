# LEGEND L200 dataflow

Implementation of an automatic data processing flow for L200
data, based on
[Snakemake](https://snakemake.readthedocs.io/).


## Configuration

Data processing resources are configured via a single site-dependent (and
possibly user-dependent) configuration file, named `config.json` in the
following. You may choose an arbitrary name, though.

Use the included [templates/config.json](templates/config.json) as a template
and adjust the data base paths as necessary. Note that, when running Snakemake,
the default path to the config file is `./config.json`.


## Key-Lists

Data generation is based on key-lists, which are flat text files
(extension ".keylist") containing one entry of the form
`{experiment}-{period}-{run}-{datatype}-{timestamp}` per line.

Key-lists can be auto-generated based on the available  DAQ files
using Snakemake targets of the form

* `all-{experiment}.keylist`
* `all-{experiment}-{period}.keylist`
* `all-{experiment}-{period}-{run}.keylist`
* `all-{experiment}-{period}-{run}-{datatype}.keylist`

which will generate the list of available file keys for all l200 files, resp.
a specific period, or a specific period and run, etc.

For example:
```shell
$ snakemake all-l200-myper.keylist
```
will generate a key-list with all files regarding period `myper`.


## File-Lists

File-lists are flat files listing output files that should be generated,
with one file per line. A file-list will typically be generated for a given
data tier from a key-list, using the Snakemake targets of the form
`{label}-{tier}.filelist` (generated from `{label}.keylist`).

For file lists based on auto-generated key-lists like
`all-{experiment}-{period}-{tier}.filelist`, the corresponding key-list
(`all-{experiment}-{period}.keylist` in this case) will be created
automatically, if it doesn't exist.

Example:
```shell
$ snakemake all-mydet-mymeas-tier2.filelist
```

File-lists may of course also be derived from custom keylists, generated
manually or by other means, e.g. `my-dataset-raw.filelist` will be
generated from `my-dataset.keylist`.


## Main output generation

Usually, the main output will be determined by a file-list, resp. a key-list
and data tier. The special output target `{label}-{tier}.gen` is used to
generate all files listed in `{label}-{tier}.filelist`. After the files
are created, the empty file `{label}-{tier}.filelist` will be created to
mark the successful data production.

Snakemake targets like `all-{experiment}-{period}-{tier}.gen` may be used
to automatically generate key-lists and file-lists (if not already present)
and produce all possible output for the given data tier, based on available
tier0 files which match the target.

Example:
```shell
$ snakemake all-mydet-mymeas-tier2.gen
```
Targets like `my-dataset-raw.gen` (derived from a key-list
`my-dataset.keylist`) are of course allowed as well.


## Monitoring

Snakemake supports monitoring by connecting to a
[panoptes](https://github.com/panoptes-organization/panoptes) server.

Run (e.g.)
```shell
$ panoptes --port 5000
```
in the background to run a panoptes server instance, which comes with a
GUI that can be accessed with a web-brower on the specified port.

Then use the Snakemake option `--wms-monitor` to instruct Snakemake to push
progress information to the panoptes server:
```shell
snakemake --wms-monitor http://127.0.0.1:5000 [...]
```

## Using software containers

This dataflow doesn't use Snakemake's internal Singularity support, but
instead supports Singularity containers via
[`venv`](https://github.com/oschulz/singularity-venv) environments
for greater control.

To use this, the path to `venv` and the name of the environment must be set
in `config.json`.

This is only relevant then running Snakemake *outside* of the software
container, e.g. then using a batch system (see below). If Snakemake
and the whole workflow is run inside of a container instance, no
container-related settings in `config.json` are required.
