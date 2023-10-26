"""
Helper functions for running data production
"""


def read_filelist(wildcards):
    with checkpoints.gen_filelist.get(
        label=wildcards.label, tier=wildcards.tier, extension="file"
    ).output[0].open() as f:
        files = f.read().splitlines()
        return files


def read_filelist_cal(wildcards, tier):
    label = f"all-{wildcards.experiment}-{wildcards.period}-{wildcards.run}-cal"
    with checkpoints.gen_filelist.get(label=label, tier=tier, extension="file").output[
        0
    ].open() as f:
        files = f.read().splitlines()
        return files


def read_filelist_pars_cal_channel(wildcards, tier):
    """
    This function will read the filelist of the channels and return a list of dsp files one for each channel
    """
    label = f"all-{wildcards.experiment}-{wildcards.period}-{wildcards.run}-cal-{wildcards.timestamp}-channels"
    with checkpoints.gen_filelist.get(label=label, tier=tier, extension="chan").output[
        0
    ].open() as f:
        files = f.read().splitlines()
        return files


def read_filelist_plts_cal_channel(wildcards, tier):
    """
    This function will read the filelist of the channels and return a list of dsp files one for each channel
    """
    label = f"all-{wildcards.experiment}-{wildcards.period}-{wildcards.run}-cal-{wildcards.timestamp}-channels"
    with checkpoints.gen_filelist.get(label=label, tier=tier, extension="chan").output[
        0
    ].open() as f:
        files = f.read().splitlines()
        files = [file.replace("par", "plt").replace("json", "pkl") for file in files]
        return files
