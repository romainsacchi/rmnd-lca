"""
This module contains constants for the filesystem paths used by Premise.
"""

from pathlib import Path

import platformdirs

# Directories for data which comes with Premise
DATA_DIR = Path(__file__).resolve().parent / "data"
INVENTORY_DIR = DATA_DIR / "additional_inventories"
VARIABLES_DIR = Path(__file__).resolve().parent / "iam_variables_mapping"
IAM_OUTPUT_DIR = DATA_DIR / "iam_output_files"

# Directories for user-created data
USER_DATA_BASE_DIR = platformdirs.user_data_path(appname="premise", appauthor="pylca")
USER_DATA_BASE_DIR.mkdir(parents=True, exist_ok=True)

DIR_CACHED_DB = USER_DATA_BASE_DIR / "cache"
DIR_CACHED_DB.mkdir(parents=True, exist_ok=True)

USER_LOGS_DIR = platformdirs.user_log_path(appname="premise", appauthor="pylca")
USER_LOGS_DIR.mkdir(parents=True, exist_ok=True)
