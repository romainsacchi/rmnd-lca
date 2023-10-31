"""
ecoinvent_modification.py exposes methods to create a database, perform transformations on it,
as well as export it back.

"""

import copy
import logging
import multiprocessing
import os
import pickle
import sys
from datetime import date
from multiprocessing import Pool as ProcessPool
from multiprocessing.pool import ThreadPool as Pool
from pathlib import Path
from typing import List, Union

import datapackage
import yaml

from . import __version__
from .cement import _update_cement
from .clean_datasets import DatabaseCleaner
from .data_collection import IAMDataCollection
from .direct_air_capture import _update_dac
from .electricity import _update_electricity
from .emissions import _update_emissions
from .export import (
    Export,
    _prepare_database,
    build_datapackage,
    generate_scenario_factor_file,
    generate_superstructure_db,
)
from .external import ExternalScenario
from .external_data_validation import check_external_scenarios, check_inventories
from .filesystem_constants import DATA_DIR, DIR_CACHED_DB, IAM_OUTPUT_DIR, INVENTORY_DIR
from .fuels import _update_fuels
from .inventory_imports import AdditionalInventory, DefaultInventory
from .report import generate_change_report, generate_summary_report
from .steel import _update_steel
from .transport import _update_vehicles
from .utils import (
    clear_existing_cache,
    create_scenario_list,
    eidb_label,
    hide_messages,
    info_on_utils_functions,
    load_constants,
    print_version,
    warning_about_biogenic_co2,
)

logger = logging.getLogger("module")

import bw2data

if int(bw2data.__version__[0]) >= 4:
    from .brightway25 import write_brightway_database

    logger.info("Using Brightway 2.5")
else:
    from .brightway2 import write_brightway_database

    logger.info("Using Brightway 2")


FILEPATH_OIL_GAS_INVENTORIES = INVENTORY_DIR / "lci-ESU-oil-and-gas.xlsx"
FILEPATH_CARMA_INVENTORIES = INVENTORY_DIR / "lci-Carma-CCS.xlsx"
FILEPATH_CO_FIRING_INVENTORIES = INVENTORY_DIR / "lci-co-firing-power-plants.xlsx"
FILEPATH_CHP_INVENTORIES = INVENTORY_DIR / "lci-combined-heat-power-plant-CCS.xlsx"
FILEPATH_CC_INVENTORIES = INVENTORY_DIR / "lci-carbon-capture.xlsx"
FILEPATH_BIOFUEL_INVENTORIES = INVENTORY_DIR / "lci-biofuels.xlsx"
FILEPATH_BIOGAS_INVENTORIES = INVENTORY_DIR / "lci-biogas.xlsx"

FILEPATH_CARBON_FIBER_INVENTORIES = INVENTORY_DIR / "lci-carbon-fiber.xlsx"
FILEPATH_HYDROGEN_DISTRI_INVENTORIES = INVENTORY_DIR / "lci-hydrogen-distribution.xlsx"

FILEPATH_HYDROGEN_INVENTORIES = INVENTORY_DIR / "lci-hydrogen-electrolysis.xlsx"
FILEPATH_HYDROGEN_SOLAR_INVENTORIES = (
    INVENTORY_DIR / "lci-hydrogen-thermochemical-water-splitting.xlsx"
)
FILEPATH_HYDROGEN_PYROLYSIS_INVENTORIES = INVENTORY_DIR / "lci-hydrogen-pyrolysis.xlsx"

FILEPATH_HYDROGEN_BIOGAS_INVENTORIES = (
    INVENTORY_DIR / "lci-hydrogen-smr-atr-biogas.xlsx"
)
FILEPATH_HYDROGEN_NATGAS_INVENTORIES = (
    INVENTORY_DIR / "lci-hydrogen-smr-atr-natgas.xlsx"
)
FILEPATH_HYDROGEN_WOODY_INVENTORIES = (
    INVENTORY_DIR / "lci-hydrogen-wood-gasification.xlsx"
)
FILEPATH_HYDROGEN_COAL_GASIFICATION_INVENTORIES = (
    INVENTORY_DIR / "lci-hydrogen-coal-gasification.xlsx"
)
FILEPATH_SYNFUEL_INVENTORIES = (
    INVENTORY_DIR / "lci-synfuels-from-FT-from-electrolysis.xlsx"
)

FILEPATH_SYNFUEL_FROM_FT_FROM_WOOD_GASIFICATION_INVENTORIES = (
    INVENTORY_DIR / "lci-synfuels-from-FT-from-wood-gasification.xlsx"
)
FILEPATH_SYNFUEL_FROM_FT_FROM_WOOD_GASIFICATION_WITH_CCS_INVENTORIES = (
    INVENTORY_DIR / "lci-synfuels-from-FT-from-wood-gasification-with-CCS.xlsx"
)
FILEPATH_SYNFUEL_FROM_FT_FROM_COAL_GASIFICATION_INVENTORIES = (
    INVENTORY_DIR / "lci-synfuels-from-FT-from-coal-gasification.xlsx"
)
FILEPATH_SYNFUEL_FROM_FT_FROM_COAL_GASIFICATION_WITH_CCS_INVENTORIES = (
    INVENTORY_DIR / "lci-synfuels-from-FT-from-coal-gasification-with-CCS.xlsx"
)

FILEPATH_SYNFUEL_FROM_BIOMASS_CCS_INVENTORIES = (
    INVENTORY_DIR / "lci-synfuels-from-FT-from-biomass-CCS.xlsx"
)
FILEPATH_SYNGAS_INVENTORIES = INVENTORY_DIR / "lci-syngas.xlsx"
FILEPATH_SYNGAS_FROM_COAL_INVENTORIES = INVENTORY_DIR / "lci-syngas-from-coal.xlsx"
FILEPATH_GEOTHERMAL_HEAT_INVENTORIES = INVENTORY_DIR / "lci-geothermal.xlsx"
FILEPATH_METHANOL_FUELS_INVENTORIES = (
    INVENTORY_DIR / "lci-synfuels-from-methanol-from-electrolysis.xlsx"
)
FILEPATH_METHANOL_CEMENT_FUELS_INVENTORIES = (
    INVENTORY_DIR / "lci-synfuels-from-methanol-from-cement-plant.xlsx"
)
FILEPATH_METHANOL_FROM_COAL_FUELS_INVENTORIES = (
    INVENTORY_DIR / "lci-synfuels-from-methanol-from-coal.xlsx"
)
FILEPATH_METHANOL_FROM_COAL_FUELS_WITH_CCS_INVENTORIES = (
    INVENTORY_DIR / "lci-synfuels-from-methanol-from-coal-with-CCS.xlsx"
)
FILEPATH_METHANOL_FROM_BIOMASS_FUELS_INVENTORIES = (
    INVENTORY_DIR / "lci-synfuels-from-methanol-from-biomass.xlsx"
)
FILEPATH_METHANOL_FROM_BIOGAS_FUELS_INVENTORIES = (
    INVENTORY_DIR / "lci-synfuels-from-methanol-from-biogas.xlsx"
)
FILEPATH_METHANOL_FROM_NATGAS_FUELS_INVENTORIES = (
    INVENTORY_DIR / "lci-synfuels-from-methanol-from-natural-gas.xlsx"
)
FILEPATH_LITHIUM = INVENTORY_DIR / "lci-lithium.xlsx"
FILEPATH_COBALT = INVENTORY_DIR / "lci-cobalt.xlsx"
FILEPATH_GRAPHITE = INVENTORY_DIR / "lci-graphite.xlsx"
FILEPATH_BATTERIES = INVENTORY_DIR / "lci-batteries.xlsx"
FILEPATH_PHOTOVOLTAICS = INVENTORY_DIR / "lci-PV.xlsx"
FILEPATH_BIGCC = INVENTORY_DIR / "lci-BIGCC.xlsx"
FILEPATH_NUCLEAR_EPR = INVENTORY_DIR / "lci-nuclear_EPR.xlsx"
FILEPATH_NUCLEAR_SMR = INVENTORY_DIR / "lci-nuclear_SMR.xlsx"
FILEPATH_WAVE = INVENTORY_DIR / "lci-wave_energy.xlsx"
FILEPATH_FUEL_CELL = INVENTORY_DIR / "lci-fuel_cell.xlsx"
FILEPATH_CSP = INVENTORY_DIR / "lci-concentrating-solar-power.xlsx"
FILEPATH_HOME_STORAGE_BATTERIES = INVENTORY_DIR / "lci-home-batteries.xlsx"
FILEPATH_VANADIUM = INVENTORY_DIR / "lci-vanadium.xlsx"
FILEPATH_VANADIUM_REDOX_BATTERY = INVENTORY_DIR / "lci-vanadium-redox-flow-battery.xlsx"
FILEPATH_HYDROGEN_TURBINE = INVENTORY_DIR / "lci-hydrogen-turbine.xlsx"

config = load_constants()


# Disable printing
def blockPrint():
    sys.stdout = open(os.devnull, "w")


# Restore printing
def enablePrint():
    sys.stdout = sys.__stdout__


def check_ei_filepath(filepath: str) -> Path:
    """Check for the existence of the file path."""

    if not Path(filepath).is_dir():
        raise FileNotFoundError(
            f"The directory for ecospold files {filepath} could not be found."
        )
    return Path(filepath)


def check_model_name(name: str) -> str:
    """Check for the validity of the IAM model name."""
    if name.lower() not in config["SUPPORTED_MODELS"]:
        raise ValueError(
            f"Only {config['SUPPORTED_MODELS']} are currently supported, not {name}."
        )
    return name.lower()


def check_pathway_name(name: str, filepath: Path, model: str) -> str:
    """Check the pathway name"""

    if name not in config["SUPPORTED_PATHWAYS"]:
        # If the pathway name is not a default one,
        # check that the filepath + pathway name
        # leads to an actual file

        if model.lower() not in name:
            name_check = "_".join((model.lower(), name))
        else:
            name_check = name

        if (filepath / name_check).with_suffix(".mif").is_file():
            return name
        if (filepath / name_check).with_suffix(".xlsx").is_file():
            return name
        if (filepath / name_check).with_suffix(".csv").is_file():
            return name
        raise ValueError(
            f"Only {config['SUPPORTED_PATHWAYS']} are currently supported, not {name}."
        )

    if model.lower() not in name:
        name_check = "_".join((model.lower(), name))
    else:
        name_check = name

    if (filepath / name_check).with_suffix(".mif").is_file():
        return name
    if (filepath / name_check).with_suffix(".xlsx").is_file():
        return name
    if (filepath / name_check).with_suffix(".csv").is_file():
        return name

    raise ValueError(
        f"Cannot find the IAM scenario file at this location: {filepath / name_check}."
    )


def check_year(year: [int, float]) -> int:
    """Check for the validity of the year passed."""
    try:
        year = int(year)
    except ValueError as err:
        raise Exception(f"{year} is not a valid year.") from err

    try:
        assert 2005 <= year <= 2100
    except AssertionError as err:
        raise Exception(f"{year} must be comprised between 2005 and 2100.") from err

    return year


def check_filepath(path: str) -> Path:
    """
    Check for the existence of the file.
    """
    if not Path(path).is_dir():
        raise FileNotFoundError(f"The filepath {path} could not be found.")
    return Path(path)


def check_exclude(list_exc: List[str]) -> List[str]:
    """
    Check for the validity of the list of excluded functions.
    """

    if not isinstance(list_exc, list):
        raise TypeError("`exclude` should be a sequence of strings.")

    if not set(list_exc).issubset(config["LIST_TRANSF_FUNC"]):
        raise ValueError(
            "One or several of the transformation that you wish to exclude is not recognized."
        )
    return list_exc


def check_additional_inventories(inventories_list: List[dict]) -> List[dict]:
    """
    Check that any additional inventories that need to be imported are properly listed.
    :param inventories_list: list of dictionaries
    :return: list of dictionaries
    """

    if not isinstance(inventories_list, list):
        raise TypeError(
            "Inventories to import need to be in a sequence of dictionaries like so:"
            "["
            "{'inventories': 'a file path', 'ecoinvent version: '3.6'},"
            " {'inventories': 'a file path', 'ecoinvent version: '3.6'}"
            "]"
        )

    for inventory in inventories_list:
        if not isinstance(inventory, dict):
            raise TypeError(
                "Inventories to import need to be in a sequence of dictionaries like so:"
                "["
                "{'inventories': 'a file path', 'ecoinvent version: '3.6'},"
                " {'inventories': 'a file path', 'ecoinvent version: '3.6'}"
                "]"
            )
        if "region_duplicate" in inventory:
            if not isinstance(inventory["region_duplicate"], bool):
                raise TypeError(
                    "`region_duplicate`must be a boolean (`True`` `False`.)"
                )

        if not all(
            i for i in inventory.keys() if i in ["inventories", "ecoinvent version"]
        ):
            raise TypeError(
                "Both `inventories` and `ecoinvent version` "
                "must be present in the list of inventories to import."
            )

        if not Path(inventory["filepath"]).is_file():
            raise FileNotFoundError(
                f"Cannot find the inventory file: {inventory['inventories']}."
            )

        if inventory["ecoinvent version"] not in config["SUPPORTED_EI_VERSIONS"]:
            raise ValueError(
                "A lot of trouble will be avoided if the additional "
                f"inventories to import are ecoinvent 3.6, 3.7, 3-8 or 3.9-compliant, not {inventory['ecoinvent version']}."
            )

    return inventories_list


def check_db_version(version: [str, float]) -> str:
    """
    Check that the ecoinvent database version is supported
    :param version:
    :return: str
    """
    version = str(version)
    if version not in config["SUPPORTED_EI_VERSIONS"]:
        raise ValueError(
            f"Only {config['SUPPORTED_EI_VERSIONS']} are currently supported, not {version}."
        )

    # convert "3.7.1" to "3.7"
    if version == "3.7.1":
        version = "3.7"

    if version == "3.9.1":
        version = "3.9"

    return version


def check_scenarios(scenario: dict, key: bytes) -> dict:
    """
    Check that the scenarios are properly formatted and that
    all the necessary info is given.
    """

    if not all(name in scenario for name in ["model", "pathway", "year"]):
        raise ValueError(
            f"Missing parameters in {scenario}. Needs to include at least `model`,"
            f"`pathway` and `year`."
        )

    if "filepath" in scenario:
        filepath = scenario["filepath"]
        scenario["filepath"] = check_filepath(filepath)
    else:
        # Note: A directory path, not a file path
        scenario["filepath"] = IAM_OUTPUT_DIR
        if key is None:
            raise ValueError(
                "You need to provide the encryption key to decrypt the IAM output files provided by `premise`."
            )

    scenario["model"] = check_model_name(scenario["model"])
    scenario["pathway"] = check_pathway_name(
        scenario["pathway"], scenario["filepath"], scenario["model"]
    )
    scenario["year"] = check_year(scenario["year"])

    if "exclude" in scenario:
        scenario["exclude"] = check_exclude(scenario["exclude"])

    return scenario


def check_system_model(system_model: str) -> str:
    """
    Check that the system model is valid.
    """

    if not isinstance(system_model, str):
        raise TypeError(
            "The argument `system_model` must be a string"
            "('consequential', 'cutoff')."
        )

    if system_model not in ("consequential", "cutoff"):
        raise ValueError(
            "The argument `system_model` must be one of the two values:"
            "'consequential', 'cutoff'."
        )

    return system_model


def check_time_horizon(time_horizon: int) -> int:
    """
    Check the validity of the time horizon provided (in years).
    :param time_horizon: time horizon (in years), to determine marginal mixes for consequential modelling.
    :return: time horizon (in years)
    """

    if time_horizon is None:
        print(
            "`time_horizon`, used to identify marginal suppliers, is not specified. "
            "It is therefore set to 20 years."
        )
        time_horizon = 20

    try:
        int(time_horizon)
    except ValueError as err:
        raise Exception(
            "`time_horizon` must be an integer with a value between 5 and 50 years."
        ) from err

    if time_horizon < 5 or time_horizon > 50:
        raise ValueError(
            "`time_horizon` must be an integer with a value between 5 and 50 years."
        )

    return int(time_horizon)


def _update_all(
    scenario,
    version,
    system_model,
    modified_datasets,
    use_absolute_efficiency,
    vehicle_type,
    gains_scenario,
):
    scenario, modified_datasets, cache = _update_vehicles(
        scenario=scenario,
        vehicle_type=vehicle_type,
        version=version,
        system_model=system_model,
        modified_datasets=modified_datasets,
    )
    scenario, modified_datasets, cache = _update_electricity(
        scenario=scenario,
        version=version,
        system_model=system_model,
        modified_datasets=modified_datasets,
        use_absolute_efficiency=use_absolute_efficiency,
        cache=cache,
    )
    scenario, modified_datasets, cache = _update_dac(
        scenario=scenario,
        version=version,
        system_model=system_model,
        modified_datasets=modified_datasets,
        cache=cache,
    )
    scenario, modified_datasets, cache = _update_cement(
        scenario=scenario,
        version=version,
        system_model=system_model,
        modified_datasets=modified_datasets,
        cache=cache,
    )
    scenario, modified_datasets, cache = _update_steel(
        scenario=scenario,
        version=version,
        system_model=system_model,
        modified_datasets=modified_datasets,
        cache=cache,
    )
    scenario, modified_datasets, cache = _update_fuels(
        scenario=scenario,
        version=version,
        system_model=system_model,
        modified_datasets=modified_datasets,
        cache=cache,
    )
    scenario, modified_datasets = _update_emissions(
        scenario, version, system_model, gains_scenario, modified_datasets
    )

    return scenario, modified_datasets


def _export_to_matrices(obj):
    obj.export_db_to_matrices()


def _export_to_simapro(obj):
    obj.export_db_to_simapro()


class NewDatabase:
    """
    Class that represents a new wurst inventory database, modified according to IAM data.

    :ivar source_type: the source of the ecoinvent database. Can be `brigthway` or `ecospold`.
    :vartype source_type: str
    :ivar source_db: name of the ecoinvent source database
    :vartype source_db: str
    :ivar source_version: version of the ecoinvent source database.
        Currently works with ecoinvent cut-off 3.5, 3.6, 3.7, 3.7.1, 3.8, 3.9 and 3.9.1.
    :vartype source_version: str
    :ivar system_model: Can be `cutoff` (default) or `consequential`.
    :vartype system_model: str

    """

    def __init__(
        self,
        scenarios: List[dict],
        source_version: str = "3.9",
        source_type: str = "brightway",
        key: bytes = None,
        source_db: str = None,
        source_file_path: str = None,
        additional_inventories: List[dict] = None,
        system_model: str = "cutoff",
        system_args: dict = None,
        use_cached_inventories: bool = True,
        use_cached_database: bool = True,
        external_scenarios: list = None,
        quiet=False,
        keep_uncertainty_data=False,
        gains_scenario="CLE",
        use_absolute_efficiency=False,
        use_multiprocessing=True,
    ) -> None:
        self.source = source_db
        self.version = check_db_version(source_version)
        self.source_type = source_type
        self.system_model = check_system_model(system_model)
        self.system_model_args = system_args
        self.use_absolute_efficiency = use_absolute_efficiency
        self.multiprocessing = use_multiprocessing

        # if version is anything other than 3.8 or 3.9
        # and system_model is "consequential"
        # raise an error
        if self.version not in ["3.8", "3.9"] and self.system_model == "consequential":
            raise ValueError(
                "Consequential system model is only available for ecoinvent 3.8 or 3.9."
            )

        if gains_scenario not in ["CLE", "MFR"]:
            raise ValueError("gains_scenario must be either 'CLE' or 'MFR'")
        self.gains_scenario = gains_scenario

        if self.source_type == "ecospold":
            self.source_file_path = check_ei_filepath(source_file_path)
        else:
            self.source_file_path = None

        self.scenarios = [check_scenarios(scenario, key) for scenario in scenarios]

        # create dictionary that keeps track of emptied and created datasets
        self.modified_datasets = {
            (s["model"], s["pathway"], s["year"]): {"emptied": [], "created": []}
            for s in self.scenarios
        }

        # print some info
        if not quiet:
            print_version()
            warning_about_biogenic_co2()
            info_on_utils_functions()
            hide_messages()

        if additional_inventories:
            self.additional_inventories = check_additional_inventories(
                additional_inventories
            )
        else:
            self.additional_inventories = None

        if external_scenarios:
            self.datapackages, self.scenarios = check_external_scenarios(
                external_scenarios, self.scenarios
            )
        else:
            self.datapackages = None

        print("\n//////////////////// EXTRACTING SOURCE DATABASE ////////////////////")
        if use_cached_database:
            self.database = self.__find_cached_db(
                source_db, keep_uncertainty_data=keep_uncertainty_data
            )
            print("Done!")
        else:
            self.database = self.__clean_database(
                keep_uncertainty_data=keep_uncertainty_data
            )

        print("\n////////////////// IMPORTING DEFAULT INVENTORIES ///////////////////")
        if use_cached_inventories:
            data = self.__find_cached_inventories(source_db)
            if data is not None:
                self.database.extend(data)

        else:
            self.__import_inventories(keep_uncertainty_data=keep_uncertainty_data)

        if self.additional_inventories:
            data = self.__import_additional_inventories(self.additional_inventories)
            self.database.extend(data)

        print("Done!")

        print("\n/////////////////////// EXTRACTING IAM DATA ////////////////////////")

        def _fetch_iam_data(scenario):
            data = IAMDataCollection(
                model=scenario["model"],
                pathway=scenario["pathway"],
                year=scenario["year"],
                external_scenarios=scenario.get("external scenarios"),
                filepath_iam_files=scenario["filepath"],
                key=key,
                system_model=self.system_model,
                system_model_args=self.system_model_args,
                gains_scenario=self.gains_scenario,
                use_absolute_efficiency=self.use_absolute_efficiency,
            )
            scenario["iam data"] = data

            if self.datapackages:
                scenario["external data"] = data.get_external_data(self.datapackages)

            scenario["database"] = copy.deepcopy(self.database)

        # use multiprocessing to speed up the process
        if self.multiprocessing:
            with Pool(processes=multiprocessing.cpu_count()) as pool:
                pool.map(_fetch_iam_data, self.scenarios)
        else:
            for scenario in self.scenarios:
                _fetch_iam_data(scenario)

        print("Done!")

    def __find_cached_db(self, db_name: str, keep_uncertainty_data: bool) -> List[dict]:
        """
        If `use_cached_db` = True, then we look for a cached database.
        If cannot be found, we create a cache for next time.
        :param db_name: database name
        :return: database
        """
        # build file path
        if db_name is None and self.source_type == "ecospold":
            db_name = f"ecospold_{self.system_model}_{self.version}"

        file_name = (
            DIR_CACHED_DB
            / f"cached_{''.join(tuple(map( str , __version__ )))}_{db_name.strip().lower()}.pickle"
        )

        # check that file path leads to an existing file
        if file_name.exists():
            # return the cached database
            return pickle.load(open(file_name, "rb"))

        # extract the database, pickle it for next time and return it
        print("Cannot find cached database. Will create one now for next time...")
        clear_existing_cache()
        database = self.__clean_database(keep_uncertainty_data=keep_uncertainty_data)
        pickle.dump(database, open(file_name, "wb"))
        return database

    def __find_cached_inventories(self, db_name: str) -> Union[None, List[dict]]:
        """
        If `use_cached_inventories` = True, then we look for a cached inventories.
        If cannot be found, we create a cache for next time.
        :param db_name: database name
        :return: database
        """
        # build file path
        if db_name is None and self.source_type == "ecospold":
            db_name = f"ecospold_{self.system_model}_{self.version}"

        file_name = (
            DIR_CACHED_DB
            / f"cached_{''.join(tuple(map( str , __version__ )))}_{db_name.strip().lower()}_inventories.pickle"
        )

        # check that file path leads to an existing file
        if file_name.exists():
            # return the cached database
            return pickle.load(open(file_name, "rb"))

        # else, extract the database, pickle it for next time and return it
        print("Cannot find cached inventories. Will create them now for next time...")
        data = self.__import_inventories()
        pickle.dump(data, open(file_name, "wb"))
        print(
            "Data cached. It is advised to restart your workflow at this point.\n"
            "This allows premise to use the cached data instead, which results in\n"
            "a faster workflow."
        )
        return None

    def __clean_database(self, keep_uncertainty_data) -> List[dict]:
        """
        Extracts the ecoinvent database, loads it into a dictionary and does a little bit of housekeeping
        (adds missing locations, reference products, etc.).
        :return:
        """
        return DatabaseCleaner(
            self.source, self.source_type, self.source_file_path, self.version
        ).prepare_datasets(keep_uncertainty_data)

    def __import_inventories(self, keep_uncertainty_data: bool = False) -> List[dict]:
        """
        This method will trigger the import of a number of pickled inventories
        and merge them into the database dictionary.
        """

        print("Importing default inventories...\n")

        # with HiddenPrints():
        # Manual import
        # file path and original ecoinvent version
        data = []
        filepaths = [
            (FILEPATH_OIL_GAS_INVENTORIES, "3.7"),
            (FILEPATH_CARMA_INVENTORIES, "3.5"),
            (FILEPATH_CO_FIRING_INVENTORIES, "3.5"),
            (FILEPATH_CHP_INVENTORIES, "3.5"),
            (FILEPATH_CC_INVENTORIES, "3.9"),
            (FILEPATH_BIOGAS_INVENTORIES, "3.6"),
            (FILEPATH_CARBON_FIBER_INVENTORIES, "3.9"),
            (FILEPATH_LITHIUM, "3.8"),
            (FILEPATH_COBALT, "3.8"),
            (FILEPATH_GRAPHITE, "3.8"),
            (FILEPATH_BATTERIES, "3.8"),
            (FILEPATH_HOME_STORAGE_BATTERIES, "3.9"),
            (FILEPATH_PHOTOVOLTAICS, "3.7"),
            (FILEPATH_HYDROGEN_INVENTORIES, "3.9"),
            (FILEPATH_HYDROGEN_SOLAR_INVENTORIES, "3.9"),
            (FILEPATH_HYDROGEN_PYROLYSIS_INVENTORIES, "3.9"),
            (FILEPATH_METHANOL_FUELS_INVENTORIES, "3.7"),
            (FILEPATH_METHANOL_CEMENT_FUELS_INVENTORIES, "3.7"),
            (FILEPATH_HYDROGEN_COAL_GASIFICATION_INVENTORIES, "3.7"),
            (FILEPATH_METHANOL_FROM_COAL_FUELS_INVENTORIES, "3.7"),
            (FILEPATH_METHANOL_FROM_COAL_FUELS_WITH_CCS_INVENTORIES, "3.7"),
            (FILEPATH_HYDROGEN_DISTRI_INVENTORIES, "3.7"),
            (FILEPATH_HYDROGEN_BIOGAS_INVENTORIES, "3.7"),
            (FILEPATH_HYDROGEN_NATGAS_INVENTORIES, "3.7"),
            (FILEPATH_HYDROGEN_WOODY_INVENTORIES, "3.7"),
            (FILEPATH_HYDROGEN_TURBINE, "3.9"),
            (FILEPATH_SYNGAS_INVENTORIES, "3.9"),
            (FILEPATH_SYNGAS_FROM_COAL_INVENTORIES, "3.7"),
            (FILEPATH_BIOFUEL_INVENTORIES, "3.7"),
            (FILEPATH_SYNFUEL_INVENTORIES, "3.7"),
            (
                FILEPATH_SYNFUEL_FROM_FT_FROM_WOOD_GASIFICATION_INVENTORIES,
                "3.7",
            ),
            (
                FILEPATH_SYNFUEL_FROM_FT_FROM_WOOD_GASIFICATION_WITH_CCS_INVENTORIES,
                "3.7",
            ),
            (
                FILEPATH_SYNFUEL_FROM_FT_FROM_COAL_GASIFICATION_INVENTORIES,
                "3.7",
            ),
            (
                FILEPATH_SYNFUEL_FROM_FT_FROM_COAL_GASIFICATION_WITH_CCS_INVENTORIES,
                "3.7",
            ),
            (FILEPATH_GEOTHERMAL_HEAT_INVENTORIES, "3.6"),
            (FILEPATH_BIGCC, "3.8"),
            (FILEPATH_NUCLEAR_EPR, "3.8"),
            (FILEPATH_NUCLEAR_SMR, "3.8"),
            (FILEPATH_WAVE, "3.8"),
            (FILEPATH_FUEL_CELL, "3.9"),
            (FILEPATH_CSP, "3.9"),
            (FILEPATH_VANADIUM, "3.8"),
            (FILEPATH_VANADIUM_REDOX_BATTERY, "3.9"),
        ]
        for filepath in filepaths:
            # make an exception for FILEPATH_OIL_GAS_INVENTORIES
            # ecoinvent version is 3.9
            if filepath[0] == FILEPATH_OIL_GAS_INVENTORIES and self.version == "3.9":
                continue

            inventory = DefaultInventory(
                database=self.database,
                version_in=filepath[1],
                version_out=self.version,
                path=filepath[0],
                system_model=self.system_model,
                keep_uncertainty_data=keep_uncertainty_data,
            )
            datasets = inventory.merge_inventory()
            data.extend(datasets)

            self.database.extend(datasets)

        # print("Done!\n")
        return data

    def __import_additional_inventories(
        self, data_package: [datapackage.DataPackage, list]
    ) -> List[dict]:
        """
        This method will trigger the import of a number of inventories
        and merge them into the database dictionary.

        :param data_package: datapackage.DataPackage or list of file paths
        :return: list of dictionaries

        """
        print("\n//////////////// IMPORTING USER-DEFINED INVENTORIES ////////////////")

        data = []

        if isinstance(data_package, list):
            # this is a list of file paths
            for file_path in data_package:
                additional = AdditionalInventory(
                    database=self.database,
                    version_in=file_path["ecoinvent version"],
                    version_out=self.version,
                    path=file_path["filepath"],
                    system_model=self.system_model,
                )
                additional.prepare_inventory()
                data.extend(additional.merge_inventory())

        elif isinstance(data_package, datapackage.DataPackage):
            if data_package.get_resource("inventories"):
                additional = AdditionalInventory(
                    database=self.database,
                    version_in=data_package.descriptor["ecoinvent"]["version"],
                    version_out=self.version,
                    path=data_package.get_resource("inventories").source,
                    system_model=self.system_model,
                )
                data.extend(additional.merge_inventory())
        else:
            raise TypeError("Unknown data type for datapackage.")

        return data

    def update_electricity(self) -> None:
        """
        This method will update the electricity inventories
        with the data from the IAM scenarios.

        """

        print("\n/////////////////////////// ELECTRICITY ////////////////////////////")

        # use multiprocessing to speed up the process
        if self.multiprocessing:
            with ProcessPool(processes=multiprocessing.cpu_count()) as pool:
                args = [
                    (
                        scenario,
                        self.version,
                        self.system_model,
                        self.modified_datasets,
                        self.use_absolute_efficiency,
                    )
                    for scenario in self.scenarios
                ]
                results = pool.starmap(_update_electricity, args)

            for s, scenario in enumerate(self.scenarios):
                self.scenarios[s] = results[s][0]
                self.modified_datasets[
                    (scenario["model"], scenario["pathway"], scenario["year"])
                ] = results[s][1][
                    (scenario["model"], scenario["pathway"], scenario["year"])
                ]
        else:
            for scenario in self.scenarios:
                scenario, self.modified_datasets, _ = _update_electricity(
                    scenario=scenario,
                    version=self.version,
                    system_model=self.system_model,
                    modified_datasets=self.modified_datasets,
                    use_absolute_efficiency=self.use_absolute_efficiency,
                )

        print("Done!\n")

    def update_dac(self) -> None:
        """
        This method will update the Direct Air Capture (DAC) inventories
        with the data from the IAM scenarios.

        """

        print("\n//////////////////////// DIRECT AIR CAPTURE ////////////////////////")

        # use multiprocessing to speed up the process
        if self.multiprocessing:
            with ProcessPool(processes=multiprocessing.cpu_count()) as pool:
                args = [
                    (
                        scenario,
                        self.version,
                        self.system_model,
                        self.modified_datasets,
                    )
                    for scenario in self.scenarios
                ]
                results = pool.starmap(_update_dac, args)

            for s, scenario in enumerate(self.scenarios):
                self.scenarios[s] = results[s][0]
                self.modified_datasets[
                    (scenario["model"], scenario["pathway"], scenario["year"])
                ] = results[s][1][
                    (scenario["model"], scenario["pathway"], scenario["year"])
                ]
        else:
            for scenario in self.scenarios:
                scenario, self.modified_datasets, _ = _update_dac(
                    scenario=scenario,
                    version=self.version,
                    system_model=self.system_model,
                    modified_datasets=self.modified_datasets,
                )

        print("Done!\n")

    def update_fuels(self) -> None:
        """
        This method will update the fuels inventories
        with the data from the IAM scenarios.
        """
        print("\n////////////////////////////// FUELS ///////////////////////////////")

        # use multiprocessing to speed up the process
        if self.multiprocessing:
            with ProcessPool(processes=multiprocessing.cpu_count()) as pool:
                args = [
                    (
                        scenario,
                        self.version,
                        self.system_model,
                        self.modified_datasets,
                    )
                    for scenario in self.scenarios
                ]
                results = pool.starmap(_update_fuels, args)

            for s, scenario in enumerate(self.scenarios):
                self.scenarios[s] = results[s][0]
                self.modified_datasets[
                    (scenario["model"], scenario["pathway"], scenario["year"])
                ] = results[s][1][
                    (scenario["model"], scenario["pathway"], scenario["year"])
                ]
        else:
            for scenario in self.scenarios:
                scenario, self.modified_datasets, _ = _update_fuels(
                    scenario=scenario,
                    version=self.version,
                    system_model=self.system_model,
                    modified_datasets=self.modified_datasets,
                )

        print("Done!\n")

    def update_cement(self) -> None:
        """
        This method will update the cement inventories
        with the data from the IAM scenarios.
        """
        print("\n///////////////////////////// CEMENT //////////////////////////////")

        # use multiprocessing to speed up the process
        if self.multiprocessing:
            with ProcessPool(processes=multiprocessing.cpu_count()) as pool:
                args = [
                    (
                        scenario,
                        self.version,
                        self.system_model,
                        self.modified_datasets,
                    )
                    for scenario in self.scenarios
                ]
                results = pool.starmap(_update_cement, args)

            for s, scenario in enumerate(self.scenarios):
                self.scenarios[s] = results[s][0]
                self.modified_datasets[
                    (scenario["model"], scenario["pathway"], scenario["year"])
                ] = results[s][1][
                    (scenario["model"], scenario["pathway"], scenario["year"])
                ]
        else:
            for scenario in self.scenarios:
                scenario, self.modified_datasets, _ = _update_cement(
                    scenario=scenario,
                    version=self.version,
                    system_model=self.system_model,
                    modified_datasets=self.modified_datasets,
                )

        print("Done!\n")

    def update_steel(self) -> None:
        """
        This method will update the steel inventories
        with the data from the IAM scenarios.
        """
        print("\n////////////////////////////// STEEL //////////////////////////////")

        # use multiprocessing to speed up the process
        if self.multiprocessing:
            with ProcessPool(processes=multiprocessing.cpu_count()) as pool:
                args = [
                    (
                        scenario,
                        self.version,
                        self.system_model,
                        self.modified_datasets,
                    )
                    for scenario in self.scenarios
                ]
                results = pool.starmap(_update_steel, args)

            for s, scenario in enumerate(self.scenarios):
                self.scenarios[s] = results[s][0]
                self.modified_datasets[
                    (scenario["model"], scenario["pathway"], scenario["year"])
                ] = results[s][1][
                    (scenario["model"], scenario["pathway"], scenario["year"])
                ]
        else:
            for scenario in self.scenarios:
                scenario, self.modified_datasets, _ = _update_steel(
                    scenario=scenario,
                    version=self.version,
                    system_model=self.system_model,
                    modified_datasets=self.modified_datasets,
                )

        print("Done!\n")

    def update_cars(self) -> None:
        """
        This method will update the cars inventories
        with the data from the IAM scenarios.
        """
        print("\n///////////////////////// PASSENGER CARS ///////////////////////////")

        # use multiprocessing to speed up the process
        if self.multiprocessing:
            with ProcessPool(processes=multiprocessing.cpu_count()) as pool:
                args = [
                    (
                        scenario,
                        "car",
                        self.version,
                        self.system_model,
                        self.modified_datasets,
                    )
                    for scenario in self.scenarios
                ]
                results = pool.starmap(_update_vehicles, args)

            for s, scenario in enumerate(self.scenarios):
                self.scenarios[s] = results[s][0]
                self.modified_datasets[
                    (scenario["model"], scenario["pathway"], scenario["year"])
                ] = results[s][1][
                    (scenario["model"], scenario["pathway"], scenario["year"])
                ]
        else:
            for scenario in self.scenarios:
                scenario, self.modified_datasets, _ = _update_vehicles(
                    scenario=scenario,
                    vehicle_type="car",
                    version=self.version,
                    system_model=self.system_model,
                    modified_datasets=self.modified_datasets,
                )

        print("Done!\n")

    def update_two_wheelers(self) -> None:
        """
        This method will update the two-wheelers inventories
        with the data from the IAM scenarios.
        """
        print("\n////////////////////////// TWO-WHEELERS ////////////////////////////")

        # use multiprocessing to speed up the process
        if self.multiprocessing:
            with ProcessPool(processes=multiprocessing.cpu_count()) as pool:
                args = [
                    (
                        scenario,
                        "two wheeler",
                        self.version,
                        self.system_model,
                        self.modified_datasets,
                    )
                    for scenario in self.scenarios
                ]
                results = pool.starmap(_update_vehicles, args)

            for s, scenario in enumerate(self.scenarios):
                self.scenarios[s] = results[s][0]
                self.modified_datasets[
                    (scenario["model"], scenario["pathway"], scenario["year"])
                ] = results[s][1][
                    (scenario["model"], scenario["pathway"], scenario["year"])
                ]
        else:
            for scenario in self.scenarios:
                scenario, self.modified_datasets, _ = _update_vehicles(
                    scenario=scenario,
                    vehicle_type="two wheeler",
                    version=self.version,
                    system_model=self.system_model,
                    modified_datasets=self.modified_datasets,
                )

        print("Done!\n")

    def update_trucks(self) -> None:
        """
        This method will update the trucks inventories
        with the data from the IAM scenarios.
        """

        print("\n////////////////// MEDIUM AND HEAVY DUTY TRUCKS ////////////////////")

        args = [
            (
                scenario,
                "truck",
                self.version,
                self.system_model,
                self.modified_datasets,
            )
            for scenario in self.scenarios
        ]

        # use multiprocessing to speed up the process
        if self.multiprocessing:
            with ProcessPool(processes=multiprocessing.cpu_count()) as pool:
                results = pool.starmap(_update_vehicles, args)

            for s, scenario in enumerate(self.scenarios):
                self.scenarios[s] = results[s][0]
                self.modified_datasets[
                    (scenario["model"], scenario["pathway"], scenario["year"])
                ] = results[s][1][
                    (scenario["model"], scenario["pathway"], scenario["year"])
                ]
        else:
            for scenario in self.scenarios:
                scenario, self.modified_datasets, _ = _update_vehicles(
                    scenario=scenario,
                    vehicle_type="truck",
                    version=self.version,
                    system_model=self.system_model,
                    modified_datasets=self.modified_datasets,
                )

        print("Done!\n")

    def update_buses(self) -> None:
        """
        This method will update the buses inventories
        with the data from the IAM scenarios.
        """

        print("\n////////////////////////////// BUSES ///////////////////////////////")

        # use multiprocessing to speed up the process
        if self.multiprocessing:
            with ProcessPool(processes=multiprocessing.cpu_count()) as pool:
                args = [
                    (
                        scenario,
                        "bus",
                        self.version,
                        self.system_model,
                        self.modified_datasets,
                    )
                    for scenario in self.scenarios
                ]
                results = pool.starmap(_update_vehicles, args)

            for s, scenario in enumerate(self.scenarios):
                self.scenarios[s] = results[s][0]
                self.modified_datasets[
                    (scenario["model"], scenario["pathway"], scenario["year"])
                ] = results[s][1][
                    (scenario["model"], scenario["pathway"], scenario["year"])
                ]
        else:
            for scenario in self.scenarios:
                scenario, self.modified_datasets, _ = _update_vehicles(
                    scenario=scenario,
                    vehicle_type="bus",
                    version=self.version,
                    system_model=self.system_model,
                    modified_datasets=self.modified_datasets,
                )

        print("Done!\n")

    def update_external_scenario(self):
        if self.datapackages:
            for i, scenario in enumerate(self.scenarios):
                if (
                    "exclude" not in scenario
                    or "update_external_scenario" not in scenario["exclude"]
                ):
                    for d, datapackage in enumerate(self.datapackages):
                        if "inventories" in [r.name for r in datapackage.resources]:
                            inventories = self.__import_additional_inventories(
                                datapackage
                            )
                        else:
                            inventories = []

                        resource = datapackage.get_resource("config")
                        config_file = yaml.safe_load(resource.raw_read())

                        checked_inventories, checked_database = check_inventories(
                            config_file,
                            inventories,
                            scenario["external data"][d],
                            scenario["database"],
                            scenario["year"],
                        )
                        scenario["database"] = checked_database
                        scenario["database"].extend(checked_inventories)

                    external_scenario = ExternalScenario(
                        database=scenario["database"],
                        model=scenario["model"],
                        pathway=scenario["pathway"],
                        iam_data=scenario["iam data"],
                        year=scenario["year"],
                        external_scenarios=self.datapackages,
                        external_scenarios_data=scenario["external data"],
                        version=self.version,
                        system_model=self.system_model,
                        modified_datasets=self.modified_datasets,
                    )
                    external_scenario.create_custom_markets()
                    scenario["database"] = external_scenario.database
            print(f"Log file of exchanges saved under {DATA_DIR / 'logs'}.")

        print("Done!\n")

    def update_emissions(self) -> None:
        """
        This method will update the hot pollutants emissions
        with the data from the GAINS model.
        """

        print("\n/////////////////////////// EMISSIONS //////////////////////////////")

        # use multiprocessing to speed up the process
        if self.multiprocessing:
            with ProcessPool(processes=multiprocessing.cpu_count()) as pool:
                args = [
                    (
                        scenario,
                        self.version,
                        self.system_model,
                        self.gains_scenario,
                        self.modified_datasets,
                    )
                    for scenario in self.scenarios
                ]
                results = pool.starmap(_update_emissions, args)

            for s, scenario in enumerate(self.scenarios):
                self.scenarios[s] = results[s][0]
                self.modified_datasets[
                    (scenario["model"], scenario["pathway"], scenario["year"])
                ] = results[s][1][
                    (scenario["model"], scenario["pathway"], scenario["year"])
                ]
        else:
            for scenario in self.scenarios:
                scenario, self.modified_datasets = _update_emissions(
                    scenario=scenario,
                    version=self.version,
                    system_model=self.system_model,
                    gains_scenario=self.gains_scenario,
                    modified_datasets=self.modified_datasets,
                )

        print("Done!\n")

    def update_all(self) -> None:
        """
        Shortcut method to execute all transformation functions.
        """

        print("`update_all()` will skip the following steps:")
        print("update_two_wheelers(), update_cars(), and update_buses()")
        print(
            "If you want to update these steps, "
            "please run them separately afterwards."
        )

        # use multiprocessing to speed up the process
        if self.multiprocessing:
            with ProcessPool(processes=multiprocessing.cpu_count()) as pool:
                args = [
                    (
                        scenario,
                        self.version,
                        self.system_model,
                        self.modified_datasets,
                        self.use_absolute_efficiency,
                        "truck",
                        self.gains_scenario,
                    )
                    for scenario in self.scenarios
                ]
                results = pool.starmap(_update_all, args)

            for s, scenario in enumerate(self.scenarios):
                self.scenarios[s] = results[s][0]
                self.modified_datasets[
                    (scenario["model"], scenario["pathway"], scenario["year"])
                ] = results[s][1][
                    (scenario["model"], scenario["pathway"], scenario["year"])
                ]
        else:
            for scenario in self.scenarios:
                scenario, self.modified_datasets = _update_all(
                    scenario=scenario,
                    version=self.version,
                    system_model=self.system_model,
                    modified_datasets=self.modified_datasets,
                    use_absolute_efficiency=self.use_absolute_efficiency,
                    vehicle_type="truck",
                    gains_scenario=self.gains_scenario,
                )

        self.update_external_scenario()

    def write_superstructure_db_to_brightway(
        self,
        name: str = f"super_db_{date.today()}",
        filepath: str = None,
        format: str = "excel",
    ) -> None:
        """
        Register a super-structure database,
        according to https://github.com/dgdekoning/brightway-superstructure
        :param name: name of the super-structure database
        :param filepath: filepath of the "scenarios difference file"
        :param format: format of the "scenarios difference file" export. Can be "excel", "csv" or "feather".
        :return: filepath of the "scenarios difference file"
        """

        if len(self.scenarios) < 2:
            raise ValueError(
                "At least two scenarios are needed to"
                "create a super-structure database."
            )

        cache = {}

        for scenario in self.scenarios:
            scenario, cache = _prepare_database(
                scenario=scenario,
                scenario_cache=cache,
                version=self.version,
                system_model=self.system_model,
                modified_datasets=self.modified_datasets,
            )

        if hasattr(self, "datapackages"):
            list_scenarios = create_scenario_list(self.scenarios, self.datapackages)
        else:
            list_scenarios = create_scenario_list(self.scenarios)

        self.database = generate_superstructure_db(
            origin_db=self.database,
            scenarios=self.scenarios,
            db_name=name,
            filepath=filepath,
            version=self.version,
            format=format,
            scenario_list=list_scenarios,
        )

        write_brightway_database(
            data=self.database,
            name=name,
        )

        # generate scenario report
        self.generate_scenario_report()
        # generate change report from logs
        self.generate_change_report()

    def write_db_to_brightway(self, name: [str, List[str]] = None):
        """
        Register the new database into an open brightway project.
        :param name: to give a (list) of custom name(s) to the database.
        Should either be a string if there's only one database to export.
        Or a list of strings if there are several databases.
        :type name: str
        """

        if name:
            if isinstance(name, str):
                name = [name]
            elif isinstance(name, list):
                if not all(isinstance(item, str) for item in name):
                    raise TypeError(
                        "`name` should be a string or a sequence of strings."
                    )
            else:
                raise TypeError("`name` should be a string or a sequence of strings.")
        else:
            name = [
                eidb_label(
                    scenario["model"],
                    scenario["pathway"],
                    scenario["year"],
                    version=self.version,
                    system_model=self.system_model,
                )
                for scenario in self.scenarios
            ]

        if len(name) != len(self.scenarios):
            raise ValueError(
                "The number of databases does not match the number of `name` given."
            )

        print("Write new database(s) to Brightway.")

        cache = {}

        for scenario in self.scenarios:
            scenario, cache = _prepare_database(
                scenario=scenario,
                scenario_cache=cache,
                version=self.version,
                system_model=self.system_model,
                modified_datasets=self.modified_datasets,
            )

        for scen, scenario in enumerate(self.scenarios):
            write_brightway_database(
                scenario["database"],
                name[scen],
            )
        # generate scenario report
        self.generate_scenario_report()
        # generate change report from logs
        self.generate_change_report()

    def write_db_to_matrices(self, filepath: str = None):
        """

        Exports the new database as a sparse matrix representation in csv files.

        :param filepath: path provided by the user to store the exported matrices.
        If it is a string, the path is used as main directory from which
        "iam model" / "pathway" / "year" subdirectories will be created.
        If it is a sequence of strings, each string becomes the directory
        under which the set of matrices is saved. If `filepath` is not provided,
        "iam model" / "pathway" / "year" subdirectories are created under
        the working directory.
        :type filepath: str or list

        """

        if filepath is not None:
            if isinstance(filepath, str):
                filepath = [
                    (Path(filepath) / s["model"] / s["pathway"] / str(s["year"]))
                    for s in self.scenarios
                ]
            elif isinstance(filepath, list):
                filepath = [Path(f) for f in filepath]
            else:
                raise TypeError(
                    f"Expected a string or a sequence of "
                    f"strings for `filepath`, not {type(filepath)}."
                )
        else:
            filepath = [
                (Path.cwd() / "export" / s["model"] / s["pathway"] / str(s["year"]))
                for s in self.scenarios
            ]

        print("Write new database(s) to matrix.")

        cache = {}

        # use multiprocessing to speed up the process
        # use multiprocessing to speed up the process

        for scenario in self.scenarios:
            scenario, cache = _prepare_database(
                scenario=scenario,
                scenario_cache=cache,
                version=self.version,
                system_model=self.system_model,
                modified_datasets=self.modified_datasets,
            )

        for scen, scenario in enumerate(self.scenarios):
            Export(scenario, filepath[scen], self.version).export_db_to_matrices()

        # generate scenario report
        self.generate_scenario_report()
        # generate change report from logs
        self.generate_change_report()

    def write_db_to_simapro(self, filepath: str = None):
        """
        Exports database as a CSV file to be imported in Simapro 9.x

        :param filepath: path provided by the user to store the exported import file
        :type filepath: str

        """

        filepath = filepath or Path(Path.cwd() / "export" / "simapro")

        if not os.path.exists(filepath):
            os.makedirs(filepath)

        print("Write Simapro import file(s).")

        cache = {}

        # use multiprocessing to speed up the process

        for scenario in self.scenarios:
            scenario, cache = _prepare_database(
                scenario=scenario,
                scenario_cache=cache,
                version=self.version,
                system_model=self.system_model,
                modified_datasets=self.modified_datasets,
            )

        for scen, scenario in enumerate(self.scenarios):
            Export(scenario, filepath, self.version).export_db_to_simapro()

        # generate scenario report
        self.generate_scenario_report()
        # generate change report from logs
        self.generate_change_report()

    def write_datapackage(self, name: str = f"datapackage_{date.today()}"):
        if not isinstance(name, str):
            raise TypeError("`name` should be a string.")

        cached_inventories = self.__find_cached_inventories(self.source)

        if not cached_inventories:
            cache_fp = DIR_CACHED_DB / f"cached_{self.source}_inventories.pickle"
            raise ValueError(f"No cached inventories found at {cache_fp}.")

        cache = {}
        # use multiprocessing to speed up the process

        for scenario in self.scenarios:
            scenario, cache = _prepare_database(
                scenario=scenario,
                scenario_cache=cache,
                version=self.version,
                system_model=self.system_model,
                modified_datasets=self.modified_datasets,
            )

        if hasattr(self, "datapackages"):
            list_scenarios = create_scenario_list(self.scenarios, self.datapackages)
        else:
            list_scenarios = create_scenario_list(self.scenarios)

        df, extra_inventories = generate_scenario_factor_file(
            origin_db=self.database,
            scenarios=self.scenarios,
            db_name=name,
            version=self.version,
            scenario_list=list_scenarios,
        )

        cached_inventories.extend(extra_inventories)

        list_scenarios = ["original"] + list_scenarios

        build_datapackage(
            df=df,
            inventories=cached_inventories,
            list_scenarios=list_scenarios,
            ei_version=self.version,
            name=name,
        )

        # generate scenario report
        self.generate_scenario_report()
        # generate change report from logs
        self.generate_change_report()

    def generate_scenario_report(
        self,
        filepath: [str, Path] = None,
        name: str = f"scenario_report_{date.today()}.xlsx",
    ):
        """
        Generate a report of the scenarios.
        """

        print("Generate scenario report.")

        if filepath is not None:
            if isinstance(filepath, str):
                filepath = Path(filepath)
        else:
            filepath = Path(Path.cwd() / "export" / "scenario_report")

        if not os.path.exists(filepath):
            os.makedirs(filepath)

        name = Path(name)
        if name.suffix != ".xlsx":
            name = name.with_suffix(".xlsx")

        generate_summary_report(self.scenarios, filepath / name)

        print(f"Report saved under {filepath}.")

    def generate_change_report(self):
        """
        Generate a report of the changes between the original database and the scenarios.
        """

        print("Generate change report.")
        generate_change_report(
            self.source, self.version, self.source_type, self.system_model
        )
        # saved under working directory
        print(f"Report saved under {os.getcwd()}.")
