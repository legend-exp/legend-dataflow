#!/usr/bin/env bash

# IMPORTANT: this script must be executed from the legend-dataflow directory

./tests/runprod/install.sh

for test in tests/runprod/test-*.sh; do
    echo "::group::test $test"
    ./"$test" || exit 1
    echo "::endgroup::"
done
