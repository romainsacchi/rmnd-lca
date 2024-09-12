"""
This module contains constants for the filesystem paths used by Premise.
"""

from pathlib import Path

import platformdirs
import yaml


def load_var_file():
    """Check if the variable file exists and load it."""
    var_file = Path.cwd() / "variables.yaml"
    if var_file.exists():
        with open(var_file, "r") as f:
            return yaml.safe_load(f)
    else:
        return None


VARIABLES = load_var_file() or {}

# Directories for data which comes with Premise
DATA_DIR = Path(__file__).resolve().parent / "data"
INVENTORY_DIR = DATA_DIR / "additional_inventories"
VARIABLES_DIR = Path(__file__).resolve().parent / "iam_variables_mapping"
IAM_OUTPUT_DIR = DATA_DIR / "iam_output_files"

if "USER_DATA_BASE_DIR" in VARIABLES:
    USER_DATA_BASE_DIR = Path(VARIABLES.get("USER_DATA_BASE_DIR"))
    print(f"USER_DATA_BASE_DIR: {USER_DATA_BASE_DIR}")
else:
    USER_DATA_BASE_DIR = platformdirs.user_data_path(
        appname="premise", appauthor="pylca"
    )
USER_DATA_BASE_DIR.mkdir(parents=True, exist_ok=True)

DIR_CACHED_DB = USER_DATA_BASE_DIR / "cache"
DIR_CACHED_DB.mkdir(parents=True, exist_ok=True)

DIR_CACHED_FILES = USER_DATA_BASE_DIR / "cached_files"
DIR_CACHED_FILES.mkdir(parents=True, exist_ok=True)

USER_LOGS_DIR = platformdirs.user_log_path(appname="premise", appauthor="pylca")
USER_LOGS_DIR.mkdir(parents=True, exist_ok=True)
