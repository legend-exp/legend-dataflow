#!/usr/bin/env bash

# IMPORTANT: this script must be *sourced* from the legend-dataflow directory

_prod_cycle="$(realpath .)"

function get_dataflow_config_value() {
    python -c "import dbetto; print(dbetto.AttrsDict(dbetto.utils.load_dict('${_prod_cycle}/dataflow-config.yaml')).${1})" \
        | sed "s|\$_|${_prod_cycle}|g"
}

run_test_command() {
    printf "\033[32m%s\033[0m\n" "INFO: running command: $*"

    output=$("$@" 2>&1)
    status=$?

    if [ $status -ne 0 ]; then
        printf "\033[31m%s\033[0m\n" "vvvvvv ERROR: command failed with status $status vvvvvv"
        echo "$output"
        printf "\033[31m%s\033[0m\n" "^^^^^^ ERROR: command failed with status $status ^^^^^^"
    fi

    return $status
}


export -f get_dataflow_config_value run_test_command

PRODENV="$(realpath ..)"
export PRODENV
