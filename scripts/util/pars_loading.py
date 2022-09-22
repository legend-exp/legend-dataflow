import snakemake as smk
import re
import os
from .utils import *
from .CalibCatalog import *
from .patterns import *
from .FileKey import *

class pars_catalog(CalibCatalog):
    
    pars_overwrite_pattern = processing_overwrite_pattern()'.{ext}'
    
    @staticmethod
    def match_pars_files(filelist1, filelist2):
        for file2 in filelist2:
            fk2 = ProcessingFileKey.get_filekey_from_pattern(file2)
            for j,file1 in enumerate(filelist1):
                fk1 = ProcessingFileKey.get_filekey_from_pattern(file1)
                if fk1.processing_step == fk2.processing_step and fk1.datatype == fk2.datatype:
                    filelist1[j]=file2
                    if len(filelist2)>1:
                        filelist2.remove(file2)
                    else:
                        filelist2 =[]
        return filelist1, filelist2
    
    @staticmethod
    def select_pars_files(filelist, processing_step):
        if isinstance(processing_step, str):
            processing_step = [processing_step]
        out_list = []
        for file in filelist:
            fk = ProcessingFileKey.get_filekey_from_pattern(file)
            if fk.processing_step in processing_step:
                out_list.append(file)
        return out_list
    
    @staticmethod
    def select_pars_overwrite_files(filelist, processing_step):
        if isinstance(processing_step, str):
            processing_step = [processing_step]
        out_list = []
        for file in filelist:
            fk = ProcessingFileKey.get_filekey_from_pattern(file, pars_catalog.pars_overwrite_pattern)
            if fk.processing_step in processing_step:
                out_list.append(file)
        return out_list
        
    @staticmethod
    def get_par_file(setup, timestamp, name):
        par_file = os.path.join(pars_path(setup),'key_resolve.jsonl')
        pars_files = pars_catalog.get_calib_files(par_file, timestamp)
        par_overwrite_file = os.path.join(par_overwrite_path(setup),'key_resolve.jsonl')
        pars_files_overwrite = pars_catalog.get_calib_files(par_overwrite_file, timestamp)
        if len(pars_files_overwrite)>0:
            pars_files, pars_files_overwrite = pars_catalog.match_pars_files(pars_files, pars_files_overwrite)
        pars_files = pars_catalog.select_pars_files(pars_files, f'par_{name}')
        pars_overwrite_files = pars_catalog.select_pars_overwrite_files(pars_files_overwrite, f'par_{name}')
        pars_files = [ProcessingFileKey.get_full_path_from_filename(par_file, f'{par_pattern()}.json', get_pattern_pars(setup, name))[0] for par_file in pars_files]
        if len(pars_overwrite_files)>0:
            pars_overwrite_files = [ProcessingFileKey.get_full_path_from_filename(par_overwrite_file, f'{par_overwrite_pattern()}.json', get_pattern_pars_overwrite(setup, name))[0] for par_overwrite_file in pars_overwrite_files]
            pars_files +=pars_overwrite_files
        return pars_files