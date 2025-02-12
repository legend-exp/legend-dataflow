#!/usr/bin/env bash

PRODENV="$(realpath ..)"
export PRODENV

python -m pip install --upgrade pip wheel setuptools
python -m pip install --upgrade '.[runprod]'

dataprod -v install --remove --system bare -- dataflow-config.yaml
