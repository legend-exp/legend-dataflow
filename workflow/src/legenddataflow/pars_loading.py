"""
This module uses the time validity resolving in calibcatalog
to determine the par and par overwrite for a particular timestamp
"""

from pathlib import Path

from dbetto.catalog import Catalog

from .FileKey import ProcessingFileKey

# from .patterns import
from .utils import get_pars_path, par_overwrite_path


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
            List of files
        filelist2 : list
            List of files

        Returns
        -------

        filelist1 : list
            List of files
        filelist2 : list
            List of files
        """
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
        return filelist1, filelist2

    @staticmethod
    def get_par_file(catalog, setup: dict, timestamp: str, tier: str) -> list:
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
        pars_files = catalog.valid_for(timestamp)
        par_overwrite_file = Path(par_overwrite_path(setup)) / tier / "validity.yaml"
        pars_files_overwrite = ParsCatalog.get_files(par_overwrite_file, timestamp)
        if len(pars_files_overwrite) > 0:
            pars_files, pars_files_overwrite = ParsCatalog.match_pars_files(
                pars_files, pars_files_overwrite
            )
        pars_files = [Path(get_pars_path(setup, tier)) / file for file in pars_files]
        if len(pars_files_overwrite) > 0:
            pars_overwrite_files = [
                Path(par_overwrite_path(setup)) / tier / file
                for file in pars_files_overwrite
            ]
            pars_files += pars_overwrite_files
        return pars_files
