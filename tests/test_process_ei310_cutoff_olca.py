import gc
import os

import bw2calc
import bw2data
import bw2io
import pytest
from dotenv import load_dotenv

from premise import NewDatabase, clear_inventory_cache
from premise.utils import delete_all_pickles

load_dotenv()

ei_user = os.environ["EI_USERNAME"]
ei_pass = os.environ["EI_PASSWORD"]
key = os.environ["IAM_FILES_KEY"]
# convert to bytes
key = key.encode()

ei_version = "3.10"
system_model = "cutoff"

scenarios = [
    {"model": "remind", "pathway": "SSP2-Base", "year": 2050},
    {"model": "image", "pathway": "SSP2-RCP19", "year": 2050},
    {"model": "tiam-ucl", "pathway": "SSP2-RCP19", "year": 2050},
]


@pytest.mark.slow
def test_brightway():
    bw2data.projects.set_current(f"ecoinvent-{ei_version}-{system_model}")
    clear_inventory_cache()

    if f"ecoinvent-{ei_version}-{system_model}" not in bw2data.databases:
        bw2io.import_ecoinvent_release(
            version=ei_version,
            system_model=system_model,
            username=ei_user,
            password=ei_pass,
        )

    ndb = NewDatabase(
        scenarios=scenarios,
        source_db=f"ecoinvent-{ei_version}-{system_model}",
        source_version=ei_version,
        key=key,
        system_model=system_model,
        biosphere_name=f"ecoinvent-{ei_version}-biosphere",
    )

    ndb.update()

    ndb.write_db_to_olca(filepath="olca_export.csv")

    # check existence of files
    assert os.path.exists("olca_export.csv")

    # destroy all objects
    del ndb
    gc.collect()
    delete_all_pickles()
