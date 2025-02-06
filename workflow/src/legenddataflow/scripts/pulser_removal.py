from pathlib import Path

import numpy as np
from dbetto.catalog import Props
from pygama.pargen.data_cleaning import get_tcm_pulser_ids


def get_pulser_mask(
    pulser_file=None, tcm_filelist=None, channel=None, pulser_multiplicity_threshold=10
):
    if pulser_file is not None:
        if not isinstance(pulser_file, list):
            pulser_file = [pulser_file]
        mask = np.array([], dtype=bool)
        for file in pulser_file:
            pulser_dict = Props.read_from(file)
            pulser_mask = np.array(pulser_dict["mask"])
            mask = np.append(mask, pulser_mask)

    elif tcm_filelist is not None:
        # get pulser mask from tcm files
        with Path(tcm_filelist).open() as f:
            tcm_files = f.read().splitlines()
        tcm_files = sorted(np.unique(tcm_files))
        _, mask = get_tcm_pulser_ids(tcm_files, channel, pulser_multiplicity_threshold)
    else:
        msg = "No pulser file or tcm filelist provided"
        raise ValueError(msg)

    return mask
