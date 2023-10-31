import os
from pathlib import Path

from premise.filesystem_constants import (
    DATA_DIR,
    DIR_CACHED_DB,
    IAM_OUTPUT_DIR,
    INVENTORY_DIR,
    VARIABLES_DIR,
)

dd_root = str(DATA_DIR.parent.parent)
cd_root = str(DIR_CACHED_DB.parent.parent)


def test_data_dir():
    assert isinstance(DATA_DIR, Path)
    assert DATA_DIR.is_dir()
    assert len(list(DATA_DIR.iterdir())) > 2
    for fp in DATA_DIR.iterdir():
        assert os.access(fp, os.R_OK)
    assert dd_root not in cd_root and cd_root not in dd_root


def test_inventory_dir():
    assert isinstance(INVENTORY_DIR, Path)
    assert INVENTORY_DIR.is_dir()
    for fp in INVENTORY_DIR.iterdir():
        assert os.access(fp, os.R_OK)
    assert len(list(INVENTORY_DIR.iterdir())) > 2
    assert str(DATA_DIR) in str(INVENTORY_DIR)


def test_variables_dir():
    assert isinstance(VARIABLES_DIR, Path)
    assert VARIABLES_DIR.is_dir()
    assert len(list(VARIABLES_DIR.iterdir())) > 2
    for fp in VARIABLES_DIR.iterdir():
        assert os.access(fp, os.R_OK)
    assert str(DATA_DIR) not in str(VARIABLES_DIR)
    assert dd_root in str(VARIABLES_DIR)
    assert cd_root not in str(VARIABLES_DIR)


def test_iam_output_dir():
    assert isinstance(IAM_OUTPUT_DIR, Path)
    assert IAM_OUTPUT_DIR.is_dir()
    for fp in IAM_OUTPUT_DIR.iterdir():
        assert os.access(fp, os.R_OK)
    assert len(list(IAM_OUTPUT_DIR.iterdir())) > 2
    assert str(DATA_DIR) in str(IAM_OUTPUT_DIR)


def test_user_data_dir():
    assert isinstance(DIR_CACHED_DB, Path)
    assert DIR_CACHED_DB.is_dir()
    assert os.access(DIR_CACHED_DB, os.W_OK)
    assert os.access(DIR_CACHED_DB, os.R_OK)
    assert cd_root in str(DIR_CACHED_DB)
    assert dd_root not in str(DIR_CACHED_DB)
    assert "cache" in str(DIR_CACHED_DB)
