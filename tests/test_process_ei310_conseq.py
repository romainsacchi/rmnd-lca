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
system_model = "consequential"


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

    if "test1" in bw2data.databases:
        del bw2data.databases["test1"]
    if "test2" in bw2data.databases:
        del bw2data.databases["test2"]
    if "test3" in bw2data.databases:
        del bw2data.databases["test3"]

    ndb.write_db_to_brightway(["test1", "test2", "test3"])

    method = [m for m in bw2data.methods if "IPCC" in m[0]][0]

    lca = bw2calc.LCA({bw2data.Database("test1").random(): 1}, method)
    lca.lci()
    lca.lcia()
    assert isinstance(lca.score, float)
    print(lca.score)

    # destroy all objects
    del ndb
    del lca
    gc.collect()
    delete_all_pickles()
