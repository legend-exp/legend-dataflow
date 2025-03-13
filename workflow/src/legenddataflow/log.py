import logging
import sys
import traceback
from logging.config import dictConfig
from pathlib import Path

from dbetto import Props


class StreamToLogger:
    """File-like stream object that redirects writes to a logger instance."""

    def __init__(self, logger, log_level=logging.ERROR):
        self.logger = logger
        self.log_level = log_level
        self.linebuf = ""

    def write(self, buf):
        for line in buf.rstrip().splitlines():
            self.logger.log(self.log_level, line.rstrip())

    def flush(self):
        pass


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

    # Redirect stderr to the logger (using the error level)
    sys.stderr = StreamToLogger(log, logging.ERROR)

    # Extract the stream from the logger's file handler.
    log_stream = None
    for handler in log.handlers:
        if hasattr(handler, "stream"):
            log_stream = handler.stream
            break
    if log_stream is None:
        log_stream = sys.stdout

    def excepthook(exc_type, exc_value, exc_traceback):
        traceback.print_exception(exc_type, exc_value, exc_traceback, file=log_stream)

    sys.excepthook = excepthook

    return log
