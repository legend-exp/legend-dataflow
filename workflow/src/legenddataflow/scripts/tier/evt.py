"""Console script ``build-tier-evt``: build the evt tier (event-level
reconstruction and cross-talk correction) from hit, dsp and tcm inputs,
optionally including ANN files."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import lh5
import numpy as np
from dbetto import AttrsDict, Props
from legenddataflowscripts.utils import (
    build_log,
    check_input_files,
    get_rule_config,
    prepare_output_paths,
)
from legendmeta import LegendMetadata
from lgdo.types import Array
from pygama.evt import build_evt

from legenddataflow.methods import ProcessingFileKey


def _resolve_channels(chmap, field, dic, timestamp):
    """Resolve one evt-config ``channels`` entry into a list of ``chXXX``
    names using the channel map valid at ``timestamp``."""
    system_map = chmap.map("system", unique=False)
    if dic["system"] not in system_map:
        msg = (
            f"channels field {field!r} requests system {dic['system']!r}, "
            f"which is not in the channel map for {timestamp}: available "
            f"systems are {sorted(system_map)}"
        )
        raise ValueError(msg)
    chans = system_map[dic["system"]]
    if "selectors" in dic:
        try:
            for k, val in dic["selectors"].items():
                chans = chans.map(k, unique=False)[val]
        except KeyError:
            chans = None
    if chans is not None:
        return [f"ch{chan}" for chan in list(chans.map("daq.rawid"))]
    return []


def build_tier_evt() -> None:
    argparser = argparse.ArgumentParser()
    argparser.add_argument("--hit-file")
    argparser.add_argument("--dsp-file")
    argparser.add_argument("--tcm-file")
    argparser.add_argument("--ann-file", nargs="*", default=[])
    argparser.add_argument("--xtc-file", nargs="*", default=[])
    argparser.add_argument("--par-files", nargs="*", default=[])

    argparser.add_argument("--datatype", required=True)
    argparser.add_argument("--timestamp", required=True)
    argparser.add_argument("--tier", required=True)

    argparser.add_argument("--configs", required=True)
    argparser.add_argument("--metadata", required=True)
    argparser.add_argument("--log")
    argparser.add_argument(
        "--buffer-len",
        type=int,
        default=None,
        help=(
            "number of events read per chunk. Overrides the 'buffer_len' option "
            "in the tier_evt rule config. Calibration files hold ~200x more "
            "events than physics files, so the per-chunk overhead dominates "
            "unless this is large enough to keep the chunk count low."
        ),
    )

    argparser.add_argument("--output")
    args = argparser.parse_args()

    if args.tier not in ("evt", "pet"):
        msg = f"unsupported tier {args.tier!r}: this script only builds 'evt' or 'pet'"
        raise ValueError(msg)

    # load in config
    df_config = get_rule_config(args.configs, "tier_evt", args.timestamp, args.datatype)
    log = build_log(df_config, args.log)

    check_input_files(args.hit_file, "--hit-file")
    check_input_files(args.dsp_file, "--dsp-file")
    check_input_files(args.tcm_file, "--tcm-file")
    check_input_files(args.par_files, "--par-files")
    check_input_files(args.xtc_file, "--xtc-file")
    check_input_files(args.ann_file, "--ann-file")
    prepare_output_paths(args.output)

    chmap = LegendMetadata(args.metadata, lazy=True).channelmap(on=args.timestamp)
    evt_config = AttrsDict(Props.read_from(df_config.inputs.evt_config))

    # CLI wins over the rule config, which wins over the pygama default
    buffer_len = args.buffer_len
    if buffer_len is None:
        buffer_len = df_config.get("options", {}).get("buffer_len", 10**4)
    buffer_len = int(buffer_len)
    log.debug("using buffer_len=%d", buffer_len)

    if args.datatype in ("phy", "xtc", "ssc", "rdc"):
        if len(args.xtc_file) == 0:
            msg = f"datatype is {args.datatype} but no xtc file was supplied"
            raise ValueError(msg)

        exp_string = evt_config.operations.geds___energy.expression
        exp_string = exp_string.replace(
            'xtalk_matrix_filename=""', f'xtalk_matrix_filename="{args.xtc_file[0]}"'
        )
        exp_string = exp_string.replace(
            'cal_par_files=""', f"cal_par_files={args.par_files}"
        )
        exp_string2 = exp_string.replace(
            'return_mode="energy"', 'return_mode="tcm_index"'
        )

        file_path_config = {
            "operations": {
                "geds___energy": {"expression": exp_string},
                "_geds___tcm_idx": {"expression": exp_string2},
            }
        }

        log.debug(json.dumps(file_path_config, indent=2))

        evt_config = Props.add_to(evt_config, file_path_config)

    # block for snakemake to fill in channel lists
    for field, dic in evt_config.channels.items():
        if isinstance(dic, dict):
            evt_config.channels[field] = _resolve_channels(
                chmap, field, dic, args.timestamp
            )

    log.debug(json.dumps(evt_config.channels, indent=2))

    evt_config["channel_mapping"] = {
        f"ch{chan}": dic.name for chan, dic in chmap.map("daq.rawid").items()
    }

    if "hardware_tcm_1" not in lh5.ls(args.tcm_file):
        msg = (
            f"'hardware_tcm_1' not found in tcm file {args.tcm_file}: "
            f"available objects are {lh5.ls(args.tcm_file)}"
        )
        raise ValueError(msg)

    file_table = {
        "tcm": (args.tcm_file, "hardware_tcm_1", "ch{}"),
        "dsp": (args.dsp_file, "dsp", "ch{}"),
        "hit": (args.hit_file, "hit", "ch{}"),
        "evt": (None, "evt"),
    }

    if len(args.ann_file) > 0:
        file_table["ann"] = (args.ann_file[0], "dsp", "ch{}")

    table = build_evt(
        file_table,
        evt_config,
        buffer_len=buffer_len,
    )

    if (
        "muon_config" in df_config.inputs
        and df_config.inputs["muon_config"] is not None
    ):
        muon_config = Props.read_from(df_config.inputs["muon_config"]["evt_config"])
        field_config = Props.read_from(df_config.inputs["muon_config"]["field_config"])
        # block for snakemake to fill in channel lists
        for field, dic in muon_config["channels"].items():
            if isinstance(dic, dict):
                muon_config["channels"][field] = _resolve_channels(
                    chmap, field, dic, args.timestamp
                )

        trigger_timestamp = table[field_config["ged_timestamp"]["table"]][
            field_config["ged_timestamp"]["field"]
        ].nda
        if "hardware_tcm_2" in lh5.ls(args.tcm_file):
            muon_table = build_evt(
                {
                    "tcm": (args.tcm_file, "hardware_tcm_2", "ch{}"),
                    "dsp": (args.dsp_file, "dsp", "ch{}"),
                    "hit": (args.hit_file, "hit", "ch{}"),
                    "evt": (None, "evt"),
                },
                muon_config,
                buffer_len=buffer_len,
            )
            lh5.write(
                obj=muon_table, name="evt_muon", lh5_file=args.output, wo_mode="a"
            )

            muon_timestamp = muon_table[field_config["muon_timestamp"]["field"]].nda
            muon_tbl_flag = muon_table[field_config["muon_flag"]["field"]].nda
            if len(muon_timestamp[muon_tbl_flag]) > 0:
                is_muon_veto_triggered = _find_matching_values_with_delay(
                    trigger_timestamp,
                    muon_timestamp[muon_tbl_flag],
                    field_config["jitter"],
                )
                muon_flag = np.isin(trigger_timestamp, is_muon_veto_triggered)
            else:
                muon_flag = np.zeros(len(trigger_timestamp), dtype=bool)
        else:
            muon_flag = np.zeros(len(trigger_timestamp), dtype=bool)
        table[field_config["output_field"]["table"]].add_column(
            field_config["output_field"]["field"], Array(muon_flag)
        )

    fk = ProcessingFileKey.get_filekey_from_pattern(Path(args.output).name)
    per = np.full(len(table), int(fk.period[1:]))
    run = np.full(len(table), int(fk.run[1:]))
    cycle = np.full(len(table), fk.timestamp, dtype="S16")
    table["trigger"].add_column("period", Array(per))
    table["trigger"].add_column("run", Array(run))
    table["trigger"].add_column("cycle", Array(cycle))

    lh5.write(obj=table, name="evt", lh5_file=args.output, wo_mode="a")


def _find_matching_values_with_delay(arr1, arr2, jit_delay):
    matching_values = []

    # Create an array with all possible delay values
    delays = np.linspace(0, int(1e9 * jit_delay), 10000) * jit_delay

    for delay in delays:
        arr2_delayed = arr2 + delay

        # Find matching values and indices
        mask = np.isin(arr1, arr2_delayed, assume_unique=True)
        matching_values.extend(arr1[mask])

    return np.unique(matching_values)
