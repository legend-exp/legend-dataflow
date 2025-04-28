"""
This module uses the time validity resolving in calibcatalog
to determine the par and par overwrite for a particular timestamp
"""

from pathlib import Path

from dbetto.catalog import Catalog

from .FileKey import FileKey, ProcessingFileKey

# from .patterns import
from .paths import det_status_path, get_pars_path, par_overwrite_path, pars_path


class ParsCatalog(Catalog):
    @staticmethod
    def match_pars_files(filelist1: list, filelist2: list) -> tuple[list, list]:
        """
        Takes 2 filelists and matches the files based on the processing step and datatype
        If the processing step and datatype are the same, the file in filelist1
        is replaced by the file in filelist2

        Parameters
        ----------

        filelist1 : list
            List of files to be overridden
        filelist2 : list
            List of files to override with

        Returns
        -------

        filelist : list
            List of files
        """
        if filelist2 is None or len(filelist2) == 0:
            return filelist1 if filelist1 is not None else []
        if filelist1 is None or len(filelist1) == 0:
            return filelist2 if filelist2 is not None else []
        for file2 in filelist2:
            fk2 = ProcessingFileKey.get_filekey_from_pattern(file2)
            for j, file1 in enumerate(filelist1):
                fk1 = ProcessingFileKey.get_filekey_from_pattern(file1)
                if (
                    fk1.processing_step == fk2.processing_step
                    and fk1.datatype == fk2.datatype
                ):
                    filelist1[j] = file2
                    if len(filelist2) > 1:
                        filelist2.remove(file2)
                    else:
                        filelist2 = []
        return filelist1

    def get_par_file(self, setup: dict, timestamp: str, tier: str) -> list:
        """
        Takes the par file and par overwrite file for a particular timestamp
        combines the two lists, applying the overwrite files to the par files
        return the list of par files

        Parameters
        ----------
        setup : dict
            Setup dictionary of paths
        timestamp : str
            Timestamp
        tier : str
            Tier of the processing step

        Returns
        -------
        list
            List of par files
        """
        allow_none = setup.get("allow_none_par", False)
        if pars_path(setup) not in get_pars_path(setup, tier):
            par_file = Path(get_pars_path(setup, tier)) / "validity.yaml"
            catalog = ParsCatalog.read_from(par_file)
            pars_files = catalog.valid_for(timestamp, allow_none=allow_none)
        else:
            pars_files = self.valid_for(timestamp, allow_none=allow_none)
        if pars_files is None:
            pars_files = []

        par_overwrite_file = Path(par_overwrite_path(setup)) / tier / "validity.yaml"
        overwrite_catalog = ParsCatalog.read_from(par_overwrite_file)
        pars_files_overwrite = overwrite_catalog.valid_for(
            timestamp, allow_none=allow_none
        )
        if pars_files_overwrite is None:
            pars_files_overwrite = []

        run_overwrite_validity = Path(det_status_path(setup)) / "run_override.yaml"
        run_overwrite_catalog = ParsCatalog.read_from(run_overwrite_validity)
        run_overwrite = run_overwrite_catalog.valid_for(timestamp, allow_none=True)
        if run_overwrite is not None and len(run_overwrite) > 0:
            run_overwrite_files = []
            for file in run_overwrite:
                fk = FileKey.get_filekey_from_pattern(file)
                run_overwrite_files.append(
                    f"{fk.datatype}/{fk.period}/{fk.run}/{file}-par_{tier}.yaml"
                )
            pars_files = ParsCatalog.match_pars_files(pars_files, run_overwrite_files)

        pars_files = [Path(get_pars_path(setup, tier)) / file for file in pars_files]
        pars_files_overwrite = [
            Path(par_overwrite_path(setup)) / tier / file
            for file in pars_files_overwrite
        ]
        pars_files += pars_files_overwrite
        return pars_files
