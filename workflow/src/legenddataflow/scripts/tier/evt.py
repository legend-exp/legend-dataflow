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


#: chunk size pygama's build_evt falls back to when nothing is configured
DEFAULT_BUFFER_LEN = 10**4


def _resolve_buffer_len(cli_value, df_config):
    """Resolve the evt chunk size: CLI over rule config over the pygama default.

    A non-integer or non-positive chunk size is rejected here: left alone it
    fails much later and far less clearly, inside the read loop.
    """
    buffer_len = cli_value
    source = "--buffer-len"
    if buffer_len is None:
        buffer_len = df_config.get("options", {}).get("buffer_len")
        source = "the tier_evt rule config"
    if buffer_len is None:
        # unset in both places, including an explicit `buffer_len: null`
        return DEFAULT_BUFFER_LEN

    try:
        buffer_len = int(buffer_len)
    except (TypeError, ValueError) as e:
        msg = f"buffer_len from {source} is not an integer: {buffer_len!r}"
        raise ValueError(msg) from e

    if buffer_len < 1:
        msg = f"buffer_len from {source} must be positive, got {buffer_len}"
        raise ValueError(msg)

    return buffer_len


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

    buffer_len = _resolve_buffer_len(args.buffer_len, df_config)
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
    """Return the values of *arr1* that follow a value of *arr2* within the window.

    The matching is one-sided by design: a germanium trigger is flagged only if
    it comes *at or after* a muon trigger, never before.

    Note that *jit_delay* is not itself the tolerance.  The window this searches
    is ``int(1e9 * jit_delay) * jit_delay`` — the jitter expressed in ns,
    truncated, times the jitter in seconds — so it grows roughly as
    ``1e9 * jit_delay**2``.  With the configured ``2.384185791015625e-07`` that
    is ~56.7 us.  The behaviour is preserved here rather than corrected, so the
    name and the quadratic scaling are worth revisiting in the config.

    The previous implementation swept 10000 delays and tested exact float
    equality with :func:`numpy.isin`.  That worked only because these unix
    timestamps and the jitter share a scale: at t ~ 1.76e9 s the float64
    spacing *is* ``2.384185791015625e-07``, so every timestamp difference is an
    integer multiple of it and the swept delays landed exactly on that grid.
    Comparing the offset directly is equivalent, does not depend on that
    coincidence, and replaces 10000 full scans with one sorted lookup.
    """
    arr1 = np.asarray(arr1)
    arr2 = np.asarray(arr2)
    if len(arr2) == 0 or len(arr1) == 0:
        return np.empty(0, dtype=arr1.dtype)

    window = int(1e9 * jit_delay) * jit_delay

    ref = np.sort(arr2)
    # the last reference value at or before each arr1 value
    idx = np.searchsorted(ref, arr1, side="right") - 1
    offset = np.where(idx >= 0, arr1 - ref[np.clip(idx, 0, None)], np.inf)

    return np.unique(arr1[(offset >= 0) & (offset <= window)])
