#!/usr/bin/env bash

# IMPORTANT: this script must be executed from the legend-dataflow directory

# shellcheck disable=SC1091

source "$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" &> /dev/null && pwd)/conftest.sh"

snakemake --workflow-profile workflow/profiles/lngs-build-raw -n all-*-daq.gen
snakemake --workflow-profile workflow/profiles/lngs-build-raw -n all-*-raw.gen
