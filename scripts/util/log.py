import logging
from logging.config import dictConfig
from pathlib import Path

from dbetto import Props


def build_log(config_dict, log_file=None):
    if "logging" in config_dict["options"]:
        log_config = config_dict["options"]["logging"]
        log_config = Props.read_from(log_config)
        if log_file is not None:
            Path(log_file).parent.mkdir(parents=True, exist_ok=True)
            log_config["handlers"]["dynamic"] = {
                "class": "logging.FileHandler",
                "level": "DEBUG",
                "formatter": "simple",
                "filename": log_file,
                "mode": "a",
            }
        dictConfig(log_config)
        log = logging.getLogger(config_dict["options"].get("logger", "prod"))
    else:
        if log_file is not None:
            Path(log_file).parent.mkdir(parents=True, exist_ok=True)
            logging.basicConfig(level=logging.INFO, filename=log_file, filemode="w")
        log = logging.getLogger(__name__)
    return log
