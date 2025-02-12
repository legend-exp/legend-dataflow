#!/usr/bin/env bash

# IMPORTANT: this script must be *sourced* from the legend-dataflow directory

_prod_cycle="$(realpath .)"

function get_dataflow_config_value() {
    python -c "import dbetto; print(dbetto.AttrsDict(dbetto.utils.load_dict('${_prod_cycle}/dataflow-config.yaml')).${1})" \
        | sed "s|\$_|${_prod_cycle}|g"
}
export -f get_dataflow_config_value

PRODENV="$(realpath ..)"
export PRODENV
