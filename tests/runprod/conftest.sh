#!/usr/bin/env bash

# IMPORTANT: this script must be *sourced* from the legend-dataflow directory

PRODCYCLE="$(realpath .)"
PRODENV="$(realpath ..)"

function get_dataflow_config_value() {
    python -c "import dbetto; print(dbetto.AttrsDict(dbetto.utils.load_dict('${PRODCYCLE}/dataflow-config.yaml')).${1})" \
        | sed "s|\$_|${PRODCYCLE}|g"
}

export PRODENV PRODCYCLE

echo "DEBUG: setting up test environment"

python -m pip --quiet install --upgrade pip wheel setuptools
python -m pip --quiet install --upgrade '.[runprod]'

dataprod -v install --remove --system bare -- dataflow-config.yaml
