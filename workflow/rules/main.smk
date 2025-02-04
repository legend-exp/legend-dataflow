import os
from datetime import datetime
from legenddataflow.utils import (
    filelist_path,
    log_path,
    tmp_par_path,
    pars_path,
    tmp_log_path,
)

timestamp = datetime.strftime(datetime.utcnow(), "%Y%m%dT%H%M%SZ")


# Create "{label}-{tier}.gen", based on "{label}.keylist" via
# "{label}-{tier}.filelist". Will implicitly trigger creation of all files in
# "{label}-{tier}.filelist". Example:
# "all[-{detector}[-{measurement}[-{run}[-{timestamp}]]]]-{tier}.gen":
rule autogen_output:
    """This is the main rule for running the data production

    The Snakemake target format is specified as follows:

        [all|sel]-{experiment}-{period}-{run}-{dataype}-{timestamp}-{tier}.gen

    where `experiment`, `period`, `run`, `datatype` can also be wildcards.
    `tier` must instead be fixed to a value. Specifying `sel` instead of `all`
    will only process data selected (valid) for analysis.

    Examples:

    > snakemake sel-l200-*-*-phy-dsp.gen
    > snakemake all-l200-p03-r001-phy-skm.gen

    At the end of the processing it will:

    - collect all warnings and errors in log files into a final summary file
    - run the `pygama.flow.FileDB` generation for output files
    - generate lists of valid keys
    """
    input:
        filelist=os.path.join(filelist_path(config), "{label}-{tier}.filelist"),
    output:
        gen_output="{label}-{tier}.gen",
        summary_log=log_path(config) + "/summary-{label}-{tier}-" + timestamp + ".log",
        warning_log=log_path(config) + "/warning-{label}-{tier}-" + timestamp + ".log",
    params:
        valid_keys_path=os.path.join(pars_path(config), "valid_keys"),
        filedb_path=os.path.join(pars_path(config), "filedb"),
        setup=lambda wildcards: config,
        basedir=basedir,
    threads: min(workflow.cores, 64)
    script:
        "../scripts/complete_run.py"
