"""
Module to create a logger with the given handler.
"""

import logging.config
from multiprocessing import Queue
from pathlib import Path

from .filesystem_constants import DATA_DIR

LOG_CONFIG = DATA_DIR / "utils" / "logging" / "logconfig.yaml"
DIR_LOG_REPORT = Path.cwd() / "export" / "logs"

# if DIR_LOG_REPORT folder does not exist
# we create it
if not Path(DIR_LOG_REPORT).exists():
    Path(DIR_LOG_REPORT).mkdir(parents=True, exist_ok=True)

# Assuming you have a global or passed-in queue for multiprocessing logging
log_queue = Queue()
is_config_loaded = False  # Flag to track if the logging config has been loaded


import logging.config
from pathlib import Path

import yaml

from .filesystem_constants import DATA_DIR

LOG_CONFIG = DATA_DIR / "utils" / "logging" / "logconfig.yaml"
DIR_LOG_REPORT = Path.cwd() / "export" / "logs"

# if DIR_LOG_REPORT folder does not exist
# we create it
if not Path(DIR_LOG_REPORT).exists():
    Path(DIR_LOG_REPORT).mkdir(parents=True, exist_ok=True)


def create_logger(handler):
    """Create a logger with the given handler."""
    with open(LOG_CONFIG, encoding="utf-8") as f:
        config = yaml.safe_load(f.read())
        logging.config.dictConfig(config)

    logger = logging.getLogger(handler)

    return logger


def empty_log_files():
    """
    Delete all log files found in DIR_LOG_REPORT.
    """
    for file in DIR_LOG_REPORT.iterdir():
        # if suffix is ".log"
        if file.suffix == ".log":
            try:
                file.unlink()
            except PermissionError:
                try:
                    # instead, let's empty the file
                    with open(file, "w") as f:
                        f.write("")
                except PermissionError:
                    pass
