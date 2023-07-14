__all__ = ("NewDatabase", "clear_cache", "get_regions_definition")
__version__ = (1, 5, 8)

from pathlib import Path

DATA_DIR = Path(__file__).resolve().parent / "data"
INVENTORY_DIR = Path(__file__).resolve().parent / "data" / "additional_inventories"
VARIABLES_DIR = Path(__file__).resolve().parent / "iam_variables_mapping"

from .ecoinvent_modification import NewDatabase
from .utils import clear_cache, get_regions_definition
