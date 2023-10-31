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
    with open(LOG_CONFIG, "r") as f:
        config = yaml.safe_load(f.read())
        logging.config.dictConfig(config)

    logger = logging.getLogger(handler)

    return logger
