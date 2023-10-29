checkpoint gen_filelist:
    """
    This rule generates the filelist. It is a checkpoint so when it is run it will update
    the dag passed on the files it finds as an output. It does this by taking in the search
    pattern, using this to find all the files that match this pattern, deriving the keys from
    the files found and generating the list of new files needed.
    """
    output:
        os.path.join(filelist_path(setup), "{label}-{tier}.{extension}list"),
    params:
        setup=lambda wildcards: setup,
        search_pattern=lambda wildcards: get_pattern(wildcards.tier),
        basedir=basedir,
        configs=configs,
        chan_maps=chan_maps,
        blinding=lambda wildcards: True if wildcards.tier == "raw" else false, 
    script:
        f"{workflow.source_path('../scripts/create_{wildcards.extension}list.py')}"


rule gen_fileDB_config:
    output:
        "fdb_config.json",
    script:
        f"{workflow.source_path('../scripts/gen_fiileDB_config.py')}"


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
        f"{workflow.source_path('../scripts/complete_run.py')} "