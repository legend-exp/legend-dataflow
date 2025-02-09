import numpy as np
from dbetto.catalog import Props


def get_pulser_mask(pulser_file):
    if not isinstance(pulser_file, list):
        pulser_file = [pulser_file]
    mask = np.array([], dtype=bool)
    for file in pulser_file:
        pulser_dict = Props.read_from(file)
        pulser_mask = np.array(pulser_dict["mask"])
        mask = np.append(mask, pulser_mask)

    return mask
