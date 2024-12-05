Developers Guide
===============

Snakemake is configured around a series of rules which specify how to generate a file/files from a set of input files.
These rules are defined in the ``Snakefile`` and in the files in the ``rules`` directory.
In general the structure is that a series of rules are defined to run on some calibration data generation
a final ``par_{tier}.yaml`` file at the end which can be used by the ``tier``` rule to generate all the files in the tier.
For most rules there are 2 versions the basic version and the partition version where the first uses a single run
while the latter will group many runs together.
This grouping is defined in the ``cal_grouping.yaml`` file in the `legend-datasets <https://github.com/legend-exp/legend-datasets>`_ repository.

Each rule has specified its inputs and outputs along with how to generate which can be
a shell command or a call to a python function. These scripts are stored in the ``scripts``` directory.
Additional parameters can also be defined.
Full details can be found at `snakemake https://snakemake.readthedocs.io/en/stable/snakefiles/rules.html)`_.
