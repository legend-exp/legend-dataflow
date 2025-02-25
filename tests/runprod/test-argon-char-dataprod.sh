#!/usr/bin/env bash

# IMPORTANT: this script must be executed from the legend-dataflow directory

# shellcheck disable=SC1091
source "$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" &> /dev/null && pwd)/conftest.sh"

rawdir="$(get_dataflow_config_value paths.tier_raw)"
mkdir -p "${rawdir}" || exit 1

function mkdir_n_touch() {
    mkdir -p "$(dirname "${1}")" || return 1
    touch "${1}" || return 1
}

rawfiles=(
    anp/p13/r002/l200-p13-r002-anp-20241217T094846Z-tier_raw.lh5
    anc/p13/r006/l200-p13-r006-anc-20241221T150249Z-tier_raw.lh5
    acs/p13/r006/l200-p13-r006-acs-20241221T150307Z-tier_raw.lh5
)

(
    cd "${rawdir}" || exit 1
    for file in "${rawfiles[@]}"; do
        mkdir_n_touch "$file"
    done
)

_smk_opts=(
    --touch
    --config allow_none_par=true
    --workflow-profile workflow/profiles/default
)

run_test_command snakemake "${_smk_opts[@]}" "all-p13-*-evt.gen" || exit 1
