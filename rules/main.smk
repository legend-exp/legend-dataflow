import os
from datetime import datetime
from scripts.util.utils import (
    filelist_path,
    log_path,
    tmp_par_path,
    pars_path,
    tmp_log_path,
)


# Create "{label}-{tier}.gen", based on "{label}.keylist" via
# "{label}-{tier}.filelist". Will implicitly trigger creation of all files
# in "{label}-{tier}.filelist".
# Example: "all[-{detector}[-{measurement}[-{run}[-{timestamp}]]]]-{tier}.gen":
rule autogen_output:
    """
    This is the main rule for running the data production,
    it is specified with:
    all-(experiment)-(period)-(run)-(dataype)-(timestamp)-'tier'.gen
    It will run the complete run script which collects all warnings
    and errors in log files into a final summary file. Also runs the file_db
    generation on new files as well as generating the json file with channels
    and fields in each file.
    """
    input:
        filelist=read_filelist,
    output:
        gen_output="{label}-{tier}.gen",
        summary_log=f"{log_path(setup)}/summary-"
        + "{label}-{tier}"
        + f"-{datetime.strftime(datetime.utcnow(), '%Y%m%dT%H%M%SZ')}.log",
        warning_log=f"{log_path(setup)}/warning-"
        + "{label}-{tier}"
        + f"-{datetime.strftime(datetime.utcnow(), '%Y%m%dT%H%M%SZ')}.log",
    params:
        log_path=tmp_log_path(setup),
        tmp_par_path=os.path.join(tmp_par_path(setup), "*_db.json"),
        valid_keys_path=os.path.join(pars_path(setup), "valid_keys"),
        filedb_path=os.path.join(pars_path(setup), "filedb"),
        setup=lambda wildcards: setup,
        basedir=basedir,
    script:
        "../scripts/complete_run.py"
