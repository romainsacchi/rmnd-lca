"""
Various utils functions.
"""

import os
import sys
from functools import lru_cache
from pathlib import Path
from typing import List

import pandas as pd
import xarray as xr
import yaml
from bw2data import databases
from bw2io.importers.base_lci import LCIImporter
from country_converter import CountryConverter
from prettytable import ALL, PrettyTable
from wurst.linking import (
    change_db_name,
    check_duplicate_codes,
    check_internal_linking,
    link_internal,
)
from wurst.searching import equals, get_many

from . import DATA_DIR, VARIABLES_DIR, __version__
from .data_collection import get_delimiter
from .geomap import Geomap

FUELS_PROPERTIES = VARIABLES_DIR / "fuels_variables.yaml"
CROPS_PROPERTIES = VARIABLES_DIR / "crops_variables.yaml"
EFFICIENCY_RATIO_SOLAR_PV = DATA_DIR / "renewables" / "efficiency_solar_PV.csv"


class HiddenPrints:
    """
    From https://stackoverflow.com/questions/8391411/how-to-block-calls-to-print
    """

    def __enter__(self):
        self._original_stdout = sys.stdout
        sys.stdout = open(os.devnull, "w", encoding="utf-8")

    def __exit__(self, exc_type, exc_val, exc_tb):
        sys.stdout.close()
        sys.stdout = self._original_stdout


def eidb_label(
    model: str, scenario: str, year: int, version: str, system_model: str = "cutoff"
) -> str:
    """
    Return a label to name a scenario.
    :param model: IAM model
    :param scenario: SSP/RCP scenario
    :param year: year
    :param version: version of the ecoinvent database
    :param system_model: cutoff or consequential
    :return: scenario label, str.
    """
    return f"ecoinvent_{system_model}_{version}_{model}_{scenario}_{year}"


def load_constants():
    """
    Load constants from the constants.yaml file.
    :return: dict
    """
    with open(VARIABLES_DIR / "constants.yaml", "r", encoding="utf-8") as stream:
        constants = yaml.safe_load(stream)

    return constants


@lru_cache
def get_fuel_properties() -> dict:
    """
    Loads a yaml file into a dictionary.
    This dictionary contains lower heating values
    CO2 emission factors, biogenic carbon share, for a number of fuel types.
    Mostly taken from ecoinvent and
    https://www.engineeringtoolbox.com/fuels-higher-calorific-values-d_169.html
    :return: dictionary that contains lower heating values
    :rtype: dict
    """

    with open(FUELS_PROPERTIES, "r", encoding="utf-8") as stream:
        fuel_props = yaml.safe_load(stream)

    return fuel_props


def get_crops_properties() -> dict:
    """
    Return a dictionary with crop names as keys and IAM labels as values
    relating to land use change CO2 per crop type
    :return: dict
    """
    with open(CROPS_PROPERTIES, "r", encoding="utf-8") as stream:
        crop_props = yaml.safe_load(stream)

    return crop_props


def get_efficiency_solar_photovoltaics() -> xr.DataArray:
    """
    Return an array with PV module efficiencies in function of year and technology.
    :return: xr.DataArray with PV module efficiencies
    """

    dataframe = pd.read_csv(
        EFFICIENCY_RATIO_SOLAR_PV, sep=get_delimiter(filepath=EFFICIENCY_RATIO_SOLAR_PV)
    )

    return dataframe.groupby(["technology", "year"]).mean()["efficiency"].to_xarray()


def default_global_location(database):
    """
    Set missing locations to ```GLO```
    for datasets in ``database``.
    Changes location if ``location``
    is missing or ``None``.
    Will add key ``location`` if missing.
    :param database: database to update
    """

    for dataset in get_many(database, *[equals("location", None)]):
        dataset["location"] = "GLO"
    return database


def get_regions_definition(model: str) -> None:
    """
    :param model: IAM model name, e.g., "remind", "image"

    Return a table containing the list of countries
    corresponding to each region of the model."""
    table = PrettyTable(["Region", "Countries"])

    geo = Geomap(model)
    country_converter = CountryConverter()

    for region in geo.iam_regions:
        list_countries = []
        for iso_2 in geo.iam_to_ecoinvent_location(region):
            if iso_2 in country_converter.ISO2["ISO2"].values:
                country_name = country_converter.convert(iso_2, to="name")
            else:
                country_name = iso_2

            list_countries.append(country_name)

        table.add_row([region, list_countries])

    table._max_width = {"Region": 50, "Countries": 125}
    table.hrules = ALL

    print(table)


def clear_existing_cache():
    """Clears the cache folder, except for files which contain __version__ in name.
    Useful when updating `premise`
    or encountering issues with
    inventories.
    """
    [
        f.unlink()
        for f in Path(DATA_DIR / "cache").glob("*")
        if f.is_file() and "__version__" not in f.name
    ]


# clear the cache folder
def clear_cache():
    [f.unlink() for f in Path(DATA_DIR / "cache").glob("*") if f.is_file()]
    print("Cache folder cleared!")


def print_version():
    print(f"premise v.{__version__}")


def info_on_utils_functions():
    """Display message to list utils functions"""

    table = PrettyTable(["Utils functions", "Description"])
    table.add_row(
        [
            "clear_cache()",
            (
                "Clears the cache folder. "
                "Useful when updating `premise`"
                "or encountering issues with "
                "inventories."
            ),
        ]
    )
    table.add_row(
        [
            "get_regions_definition(model)",
            "Retrieves the list of countries for each region of the model.",
        ]
    )
    table.add_row(
        [
            "ndb.NewDatabase(...)\nndb.generate_scenario_report()",
            "Generates a summary of the most important scenarios' variables.",
        ]
    )
    # align text to the left
    table.align = "l"
    table.hrules = ALL
    table._max_width = {"Utils functions": 50, "Description": 32}
    print(table)


def check_database_name(data: List[dict], name: str) -> List[dict]:
    for ds in data:
        ds["database"] = name

        for exc in ds["exchanges"]:
            if exc["type"] in ["production", "technosphere"]:
                if "input" in exc:
                    del exc["input"]

    return data


def warning_about_biogenic_co2() -> None:
    """
    Prints a simple warning about characterizing biogenic CO2 flows.
    :return: Does not return anything.
    """
    table = PrettyTable(["Warning"])
    table.add_row(
        [
            "Because some of the scenarios can yield LCI databases\n"
            "containing net negative emission technologies (NET),\n"
            "it is advised to account for biogenic CO2 flows when calculating\n"
            "Global Warming potential indicators.\n"
            "`premise_gwp` provides characterization factors for such flows.\n"
            "It also provides factors for hydrogen emissions to air.\n\n"
            "Within your bw2 project:\n"
            "from premise_gwp import add_premise_gwp\n"
            "add_premise_gwp()"
        ]
    )
    # align text to the left
    table.align = "l"
    print(table)


def hide_messages():
    """
    Hide messages from the console.
    """

    print("Keep uncertainty data?")
    print("NewDatabase(..., keep_uncertainty_data=True)")
    print("")
    print("Hide these messages?")
    print("NewDatabase(..., quiet=True)")


class PremiseImporter(LCIImporter):
    def __init__(self, db_name, data):
        self.db_name = db_name
        self.data = data
        for act in self.data:
            act["database"] = self.db_name

    # we override `write_database`
    # to allow existing databases
    # to be overwritten
    def write_database(self):
        if self.db_name in databases:
            print(f"Database {self.db_name} already exists: " "it will be overwritten.")
        super().write_database()


def clean_up(exc):
    """Remove keys from ``exc`` that are not in the schema."""

    FORBIDDEN_FIELDS_TECH = [
        "categories",
    ]

    FORBIDDEN_FIELDS_BIO = ["location", "product"]

    for field in list(exc.keys()):
        if exc[field] is None or exc[field] == "None":
            del exc[field]
            continue

        if exc["type"] == "biosphere" and field in FORBIDDEN_FIELDS_BIO:
            del exc[field]
        if exc["type"] == "technosphere" and field in FORBIDDEN_FIELDS_TECH:
            del exc[field]

    return exc


def write_brightway2_database(data, name):
    # Restore parameters to Brightway2 format
    # which allows for uncertainty and comments

    for ds in data:
        if "parameters" in ds:
            if not isinstance(ds["parameters"], list):
                if isinstance(ds["parameters"], dict):
                    ds["parameters"] = [
                        {"name": k, "amount": v} for k, v in ds["parameters"].items()
                    ]
                else:
                    ds["parameters"] = [ds["parameters"]]
            else:
                ds["parameters"] = [
                    {"name": k, "amount": v}
                    for o in ds["parameters"]
                    for k, v in o.items()
                ]

        for key, value in list(ds.items()):
            if not value:
                del ds[key]

        ds["exchanges"] = [clean_up(exc) for exc in ds["exchanges"]]

    change_db_name(data, name)
    link_internal(data)
    check_internal_linking(data)
    check_duplicate_codes(data)
    PremiseImporter(name, data).write_database()


def delete_log():
    """
    Delete log file.
    It is located in the working directory.
    """
    log_path = Path.cwd() / "premise.log"
    if log_path.exists():
        log_path.unlink()
