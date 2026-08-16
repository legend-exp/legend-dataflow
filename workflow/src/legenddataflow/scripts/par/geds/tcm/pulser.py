"""Console script ``par-geds-tcm-pulser``: identify pulser events in the tcm
tier and save their ids."""

from __future__ import annotations

import argparse

from dbetto.catalog import Props
from legenddataflowscripts.utils import (
    build_log,
    expand_filelist,
    get_rule_config,
    prepare_output_paths,
)
from pygama.pargen.data_cleaning import get_tcm_pulser_ids


def par_geds_tcm_pulser() -> None:
    argparser = argparse.ArgumentParser()
    argparser.add_argument("--configs", help="configs path", type=str, required=True)
    argparser.add_argument("--metadata", help="metadata", type=str, required=True)
    argparser.add_argument("--log", help="log file", type=str)

    argparser.add_argument("--datatype", help="Datatype", type=str, required=True)
    argparser.add_argument("--timestamp", help="Timestamp", type=str, required=True)
    argparser.add_argument("--channel", help="Channel", type=str, required=True)
    argparser.add_argument("--rawid", help="rawid", type=str, required=True)

    argparser.add_argument("--pulser-file", help="pulser file", type=str, required=True)

    argparser.add_argument("--tcm-files", help="tcm_files", nargs="*", type=str)
    args = argparser.parse_args()

    config_dict = get_rule_config(
        args.configs, "pars_tcm_pulser", args.timestamp, args.datatype
    )

    build_log(config_dict, args.log)

    prepare_output_paths(args.pulser_file)

    kwarg_dict = config_dict["inputs"]["pulser_config"]
    kwarg_dict = Props.read_from(kwarg_dict)

    # get pulser mask from tcm files
    tcm_files = expand_filelist(args.tcm_files, "--tcm-files")
    ids, mask = get_tcm_pulser_ids(
        tcm_files, args.rawid, kwarg_dict.pop("pulser_multiplicity_threshold")
    )

    Props.write_to(args.pulser_file, {"idxs": ids.tolist(), "mask": mask.tolist()})
