"""
This module uses the time validity resolving in calibcatalog
to determine the par and par overwrite for a particular timestamp
"""

import os

from .catalog import Catalog
from .FileKey import ProcessingFileKey

# from .patterns import
from .utils import get_pars_path, par_overwrite_path


class pars_catalog(Catalog):
    @staticmethod
    def match_pars_files(filelist1, filelist2):
        for file2 in filelist2:
            fk2 = ProcessingFileKey.get_filekey_from_pattern(file2)
            for j, file1 in enumerate(filelist1):
                fk1 = ProcessingFileKey.get_filekey_from_pattern(file1)
                if fk1.processing_step == fk2.processing_step and fk1.datatype == fk2.datatype:
                    filelist1[j] = file2
                    if len(filelist2) > 1:
                        filelist2.remove(file2)
                    else:
                        filelist2 = []
        return filelist1, filelist2

    @staticmethod
    def get_par_file(setup, timestamp, tier):
        par_file = os.path.join(get_pars_path(setup, tier), "validity.yaml")
        pars_files = pars_catalog.get_calib_files(par_file, timestamp)
        par_overwrite_file = os.path.join(par_overwrite_path(setup), tier, "validity.yaml")
        pars_files_overwrite = pars_catalog.get_calib_files(par_overwrite_file, timestamp)
        if len(pars_files_overwrite) > 0:
            pars_files, pars_files_overwrite = pars_catalog.match_pars_files(
                pars_files, pars_files_overwrite
            )
        pars_files = [os.path.join(get_pars_path(setup, tier), file) for file in pars_files]
        if len(pars_files_overwrite) > 0:
            pars_overwrite_files = [
                os.path.join(par_overwrite_path(setup), tier, file)
                for file in pars_files_overwrite
            ]
            pars_files += pars_overwrite_files
        return pars_files
