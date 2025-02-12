#!/usr/bin/env bash

# IMPORTANT: this script must be executed from the legend-dataflow directory

echo "DEBUG: setting up test environment"

PRODENV="$(realpath ..)"
export PRODENV

python -m pip --quiet install --upgrade pip wheel setuptools
python -m pip --quiet install --upgrade '.[runprod]'

dataprod -v install --remove --system bare -- dataflow-config.yaml
