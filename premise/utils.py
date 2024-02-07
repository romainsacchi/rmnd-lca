"""
Various utils functions.
"""

import os
import sys
import uuid
from functools import lru_cache
from numbers import Number
from pathlib import Path
from typing import Optional

import pandas as pd
import xarray as xr
import yaml
from country_converter import CountryConverter
from prettytable import ALL, PrettyTable
from wurst import rescale_exchange
from wurst.searching import biosphere, equals, get_many, technosphere

from . import __version__
from .data_collection import get_delimiter
from .filesystem_constants import DATA_DIR, DIR_CACHED_DB, VARIABLES_DIR
from .geomap import Geomap

FUELS_PROPERTIES = VARIABLES_DIR / "fuels_variables.yaml"
CROPS_PROPERTIES = VARIABLES_DIR / "crops_variables.yaml"
EFFICIENCY_RATIO_SOLAR_PV = DATA_DIR / "renewables" / "efficiency_solar_PV.csv"


def rescale_exchanges(
    ds,
    value,
    technosphere_filters=None,
    biosphere_filters=None,
    remove_uncertainty=False,
):
    """
    Adapted from wurst's change_exchanges_by_constant_factor
    but adds the possibility to preserve uncertainty data.
    """
    assert isinstance(ds, dict), "Must pass dataset dictionary document"
    assert isinstance(value, Number), "Constant factor ``value`` must be a number"

    for exc in technosphere(ds, *(technosphere_filters or [])):
        rescale_exchange(exc, value, remove_uncertainty)
    for exc in biosphere(ds, *(biosphere_filters or [])):
        rescale_exchange(exc, value, remove_uncertainty)

    return ds


# Disable printing
def blockPrint():
    with open(os.devnull, "w") as devnull:
        sys.stdout = devnull


# Restore printing
def enablePrint():
    sys.stdout = sys.__stdout__


class HiddenPrints:
    """
    From https://stackoverflow.com/questions/8391411/how-to-block-calls-to-print
    """

    def __init__(self):
        self._original_stdout = None

    def __enter__(self):
        self._original_stdout = sys.stdout
        sys.stdout = open(os.devnull, "w", encoding="utf-8")

    def __exit__(self, exc_type, exc_val, exc_tb):
        sys.stdout.close()
        sys.stdout = self._original_stdout


def eidb_label(
    model: str,
    scenario: str,
    year: int,
    version: str,
    system_model: str = "cutoff",
    datapackages=None,
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

    if datapackages is None:
        name = f"ecoinvent_{system_model}_{version}_{model}_{scenario}_{year}"
    else:
        name = f"ecoinvent_{system_model}_{version}_{model}_{scenario}_{year}_{'_'.join([dp.descriptor.get('name', 'External model') for dp in datapackages])}"
    return name


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


def clear_existing_cache(all_versions: Optional[bool] = False) -> None:
    """Clears the cache folder, except for files which contain __version__ in name.
    Useful when updating `premise`
    or encountering issues with
    inventories.
    """
    [
        f.unlink()
        for f in DIR_CACHED_DB.glob("*")
        if f.is_file()
        and (all_versions or "".join(tuple(map(str, __version__))) not in f.name)
    ]


# clear the cache folder
def clear_cache() -> None:
    clear_existing_cache(all_versions=True)
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
    print("Disable multiprocessing?")
    print("NewDatabase(..., use_multiprocessing=False)")
    print("")
    print("Hide these messages?")
    print("NewDatabase(..., quiet=True)")


def reset_all_codes(data):
    """
    Re-generate all codes in each dataset of a database
    Remove all code for each production and technosphere exchanges
    in each dataset.
    """
    for ds in data:
        ds["code"] = str(uuid.uuid4())
        for exc in ds["exchanges"]:
            if exc["type"] in ["production", "technosphere"]:
                if "input" in exc:
                    del exc["input"]

    return data


def delete_log():
    """
    Delete log file.
    It is located in the working directory.
    """
    log_path = Path.cwd() / "premise.log"
    if log_path.exists():
        log_path.unlink()


def create_scenario_list(scenarios, datapackages=None):
    list_scenarios = [f"{s['model']} - {s['pathway']} - {s['year']}" for s in scenarios]

    if "external scenarios" in scenarios[0]:
        external_model_name = "External model"
        for s, scenario in enumerate(scenarios):
            for e, ext_scenario in enumerate(scenario["external scenarios"]):
                if datapackages is not None:
                    external_model_name = datapackages[e].descriptor.get(
                        "name", "External model"
                    )
                list_scenarios[s] += f" - {external_model_name} - {ext_scenario}"

    return list_scenarios
