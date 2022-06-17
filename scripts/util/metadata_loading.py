from .Props import *
from .CalibCatalog import *
from .patterns import *


class config_catalog(CalibCatalog):
    @staticmethod
    def get_config(jsonl_file, config_path, timestamp, category):
    
        channel_path = config_catalog.get_calib_files(jsonl_file, timestamp, category)

        if isinstance(channel_path, list):
            channel_files = [os.path.join(config_path,chan) for chan in channel_path]
        else:
            channel_files = [os.path.join(config_path,channel_path)]
        channel_dict = Props.read_from(channel_files, subst_pathvar = True)
        return channel_dict
