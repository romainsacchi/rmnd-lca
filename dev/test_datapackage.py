import bw2data
from datapackage import Package

from premise import *

bw2data.projects.set_current("ei39")
sps = Package("/Users/romain/GitHub/sweet_sure-2050-switzerland/datapackage.json")
ndb = NewDatabase(
    scenarios=[
        {"model": "image", "pathway": "SSP2-RCP26", "external scenarios": [{"scenario": "SPS1", "data": sps}], "year": 2020},
    ],
    source_db="ecoinvent 3.9.1 cutoff", # <-- name of the database in the BW2 project. Must be a string.
    source_version="3.9", # <-- version of ecoinvent. Can be "3.5", "3.6", "3.7" or "3.8". Must be a string.
    key="tUePmX_S5B8ieZkkM7WUU2CnO8SmShwmAeWK9x2rTFo=",
    use_multiprocessing=False
)

ndb.update()