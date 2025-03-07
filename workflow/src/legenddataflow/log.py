import logging
from logging.config import dictConfig
from pathlib import Path

from dbetto import Props


def build_log(
    config_dict: dict, log_file: str | None = None, fallback: str = "prod"
) -> logging.Logger:
    """Build a logger from a configuration dictionary.

    If a log file is provided, the logger will write to that file.

    Parameters
    ----------
    config_dict
        A dictionary containing the logging configuration.
    log_file
        The path to the log file.
    """
    if "logging" in config_dict["options"]:
        log_config = config_dict["options"]["logging"]
        log_config = Props.read_from(log_config)

        if log_file is not None:
            Path(log_file).parent.mkdir(parents=True, exist_ok=True)
            log_config["handlers"]["dataflow"]["filename"] = log_file

        dictConfig(log_config)
        log = logging.getLogger(config_dict["options"].get("logger", "prod"))

    else:
        if log_file is not None:
            Path(log_file).parent.mkdir(parents=True, exist_ok=True)
            logging.basicConfig(level=logging.INFO, filename=log_file, filemode="w")

        log = logging.getLogger(fallback)

    return log
