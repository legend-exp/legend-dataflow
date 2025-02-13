#!/usr/bin/env bash

# IMPORTANT: this script must be executed from the legend-dataflow directory

./tests/runprod/install.sh

for test in tests/runprod/test-*.sh; do
    printf "\033[32m%s\033[0m\n" "INFO: running test $test"
    ./"$test"
done
