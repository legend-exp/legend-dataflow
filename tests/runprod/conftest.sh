#!/usr/bin/env bash

PRODENV="$(realpath ..)"
export PRODENV

echo "INFO: setting up test environment"

python -m pip --quiet install --upgrade pip wheel setuptools
python -m pip --quiet install --upgrade '.[runprod]'

dataprod -v install --remove --system bare -- dataflow-config.yaml
