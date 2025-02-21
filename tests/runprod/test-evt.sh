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
    phy/p04/r001/l200-p04-r001-phy-20230421T174901Z-tier_raw.lh5
    phy/p04/r000/l200-p04-r000-phy-20230415T033517Z-tier_raw.lh5
    phy/p03/r001/l200-p03-r001-phy-20230318T015140Z-tier_raw.lh5
    phy/p03/r000/l200-p03-r000-phy-20230312T043356Z-tier_raw.lh5
    phy/p03/r002/l200-p03-r002-phy-20230324T205907Z-tier_raw.lh5
    cal/p04/r001/l200-p04-r001-cal-20230421T131817Z-tier_raw.lh5
    cal/p04/r000/l200-p04-r000-cal-20230414T215158Z-tier_raw.lh5
    cal/p03/r001/l200-p03-r001-cal-20230317T211819Z-tier_raw.lh5
    cal/p03/r000/l200-p03-r000-cal-20230311T235840Z-tier_raw.lh5
    cal/p03/r002/l200-p03-r002-cal-20230324T161401Z-tier_raw.lh5
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

inputs="$(get_dataflow_config_value paths.metadata)"

# FIXME: remove these at some point
touch "$inputs/dataprod/overrides/dsp/cal/p03/r000/l200-p03-r000-cal-20230311T235840Z-par_dsp_svm_train.lh5"
touch "$inputs/dataprod/overrides/dsp/cal/p04/r000/l200-p04-r000-cal-20230414T215158Z-par_dsp_svm_train.lh5"

_smk_opts=(
    --touch
    --workflow-profile workflow/profiles/default
)

run_test_command snakemake "${_smk_opts[@]}" "all-*-evt.gen" || exit 1
