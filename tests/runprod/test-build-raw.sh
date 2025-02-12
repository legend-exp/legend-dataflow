#!/usr/bin/env bash

# IMPORTANT: this script must be executed from the legend-dataflow directory

# shellcheck disable=SC1091
source "$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" &> /dev/null && pwd)/conftest.sh"

sandbox=$(get_dataflow_config_value paths.sandbox_path)
mkdir -p "${sandbox}"

(
    cd "${sandbox}" || exit 1
    touch \
        l200-p03-r000-cal-20230311T235840Z.orca \
        l200-p03-r001-cal-20230317T211819Z.orca \
        l200-p03-r002-cal-20230324T161401Z.orca \
        l200-p04-r000-cal-20230414T215158Z.orca \
        l200-p04-r001-cal-20230421T131817Z.orca \
        l200-p03-r000-phy-20230312T043356Z.orca \
        l200-p03-r001-phy-20230318T015140Z.orca \
        l200-p03-r002-phy-20230324T205907Z.orca \
        l200-p04-r000-phy-20230415T033517Z.orca \
        l200-p04-r001-phy-20230421T174901Z.orca \
        l200-p13-r006-acs-20241221T150307Z.fcio \
        l200-p13-r006-anc-20241221T150249Z.fcio \
        l200-p13-r002-anp-20241217T094846Z.fcio
)

_smk_opts=(
    --workflow-profile workflow/profiles/lngs-build-raw
    --config system=bare
    --touch
)

for tier in daq raw; do
    snakemake "${_smk_opts[@]}" "all-*-${tier}.gen"
done

rm -rf "${sandbox}"
