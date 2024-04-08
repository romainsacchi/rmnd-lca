"""
data_collection.py contains the IAMDataCollection class which collects a number of data,
mostly from the IAM file. This class will have offer market shares, efficiency
and emission values for different sectors, carbon capture rates, etc.
"""

import copy
import csv
import os
from functools import lru_cache
from io import BytesIO, StringIO
from itertools import chain
from pathlib import Path
from typing import Dict, List, Union

import numpy as np
import pandas as pd
import xarray as xr
import yaml
from cryptography.fernet import Fernet
from prettytable import PrettyTable

from .filesystem_constants import DATA_DIR, IAM_OUTPUT_DIR, VARIABLES_DIR
from .geomap import Geomap
from .marginal_mixes import consequential_method

IAM_ELEC_VARS = VARIABLES_DIR / "electricity_variables.yaml"
IAM_FUELS_VARS = VARIABLES_DIR / "fuels_variables.yaml"
IAM_BIOMASS_VARS = VARIABLES_DIR / "biomass_variables.yaml"
IAM_CROPS_VARS = VARIABLES_DIR / "crops_variables.yaml"
IAM_CEMENT_VARS = VARIABLES_DIR / "cement_variables.yaml"
IAM_STEEL_VARS = VARIABLES_DIR / "steel_variables.yaml"
IAM_DAC_VARS = VARIABLES_DIR / "direct_air_capture_variables.yaml"
IAM_OTHER_VARS = VARIABLES_DIR / "other_variables.yaml"
FILEPATH_FLEET_COMP = IAM_OUTPUT_DIR / "fleet_files" / "fleet_all_vehicles.csv"
FILEPATH_IMAGE_TRUCKS_FLEET_COMP = (
    IAM_OUTPUT_DIR / "fleet_files" / "image_fleet_trucks.csv"
)
VEHICLES_MAP = DATA_DIR / "transport" / "vehicles_map.yaml"
IAM_CARBON_CAPTURE_VARS = VARIABLES_DIR / "carbon_capture_variables.yaml"
CROPS_PROPERTIES = VARIABLES_DIR / "crops_variables.yaml"
GAINS_GEO_MAP = VARIABLES_DIR / "gains_regions_mapping.yaml"
COAL_POWER_PLANTS_DATA = DATA_DIR / "electricity" / "coal_power_emissions_2012_v1.csv"


def print_missing_variables(missing_vars):
    if missing_vars:
        print("The following variables are missing from the IAM file:")
    table = PrettyTable(
        [
            "Variable",
        ]
    )
    for v in missing_vars:
        table.add_row([v])
    print(table)


def get_delimiter(data=None, filepath=None) -> str:
    sniffer = csv.Sniffer()
    if filepath:
        with open(filepath, "r", encoding="utf-8") as stream:
            data = stream.readline()
    delimiter = str(sniffer.sniff(data).delimiter)
    return delimiter


def get_crops_properties() -> dict:
    """
    Return a dictionary with crop names as keys and IAM labels as values
    relating to land use change CO2 per crop type
    :return: dict
    """
    with open(CROPS_PROPERTIES, "r", encoding="utf-8") as stream:
        crop_props = yaml.safe_load(stream)

    return crop_props


def get_oil_product_volumes(model) -> pd.DataFrame:
    """
    Load the file `oil_product_volumes.csv` that contains recent oil product volumes
    consumed in each country. This file is used to estimate the split of gasoline, diesel, LPG and kerosene
    from the `liquid fossil` category in the IAM data.
    """
    filepath = DATA_DIR / "fuels" / "oil_product_volumes.csv"
    df = pd.read_csv(
        filepath,
        delimiter=get_delimiter(filepath=filepath),
        low_memory=False,
        na_filter=False,
    )

    geo = Geomap(model=model)

    df["region"] = df.apply(lambda x: geo.ecoinvent_to_iam_location(x.country), axis=1)
    df = df.drop("country", axis=1)
    df = df.groupby(["region"]).sum()
    # add the world
    df.loc["World"] = df.sum()
    # normalize by the sum
    df = df.div(df.sum(axis=1), axis=0)

    return df


def get_metals_intensity_factors_data() -> xr.DataArray:
    """
    Read the materials intensity factors csv file and return an `xarray` with dimensions:

    * year
    * metal
    * origin_val
    * variable

    This data is further used in metals.py.
    """

    filepath = Path(DATA_DIR / "metals" / "metals_db.csv")
    df = pd.read_csv(filepath)
    df = df.melt(
        id_vars=["metal", "year", "origin_var"],
        value_vars=["mean", "median", "min", "max"],
    )

    array = (
        df.groupby(["metal", "origin_var", "year", "variable"])
        .mean()["value"]
        .to_xarray()
    )

    array = array.interpolate_na(dim="year", method="nearest", fill_value="extrapolate")
    array = array.bfill(dim="year")
    array = array.ffill(dim="year")
    array = array.fillna(0)

    return array


@lru_cache
def get_gains_IAM_data(model, gains_scenario):
    filepath = Path(
        DATA_DIR / "GAINS_emission_factors" / "iam_data" / gains_scenario
    ).glob("*")

    list_arrays = []

    for file in filepath:
        df = pd.read_csv(
            file, sep=get_delimiter(filepath=file), encoding="utf-8", low_memory=False
        )
        df = df.rename(columns={"Region": "region", "EMF30 Sector": "sector"})
        df = df.rename(columns={str(v): int(v) for v in range(1990, 2055, 5)})

        df["pollutant"] = file.stem

        array = (
            df.melt(
                id_vars=["region", "sector", "pollutant"],
                value_vars=range(1990, 2055, 5),
                var_name="year",
                value_name="value",
            )
            .groupby(["region", "sector", "year", "pollutant"])["value"]
            .mean()
            .to_xarray()
        )

        array = array.interpolate_na(
            dim="year", method="nearest", fill_value="extrapolate"
        )
        array = array.bfill(dim="year")
        array = array.ffill(dim="year")

        list_arrays.append(array)

    arr = xr.concat(list_arrays, dim="pollutant")

    with open(GAINS_GEO_MAP, "r", encoding="utf-8") as stream:
        geo_map = yaml.safe_load(stream)

    arr.coords["region"] = [geo_map[v][model] for v in arr.region.values]
    arr = arr.drop_duplicates(dim="region")

    return arr


@lru_cache
def get_gains_EU_data() -> xr.DataArray:
    """
    Read the GAINS emissions csv file and return an `xarray` with dimensions:

    * region
    * pollutant
    * sector
    * year

    :return: a multidimensional array with GAINS emissions data

    """

    filename = "GAINS_emission_factors_EU.csv"
    filepath = DATA_DIR / "GAINS_emission_factors" / filename

    gains_emi_EU = pd.read_csv(
        filepath,
        delimiter=get_delimiter(filepath=filepath),
        low_memory=False,
        dtype={
            "Region": str,
            "Sector": str,
            "Activity": str,
            "variable": str,
            "value": float,
            "year": int,
            "substance": str,
            "Activity_long": str,
        },
        encoding="utf-8",
    )
    gains_emi_EU["sector"] = gains_emi_EU["Sector"] + gains_emi_EU["Activity"]
    gains_emi_EU.drop(
        [
            "Sector",
            "Activity",
        ],
        axis=1,
    )

    gains_emi_EU = gains_emi_EU[~gains_emi_EU["value"].isna()]

    gains_emi_EU = gains_emi_EU.rename(
        columns={"Region": "region", "substance": "pollutant"}
    )

    array = (
        gains_emi_EU.groupby(["region", "pollutant", "year", "sector"])["value"]
        .mean()
        .to_xarray()
    )

    array = array.interpolate_na(dim="year", method="nearest", fill_value="extrapolate")
    array = array.bfill(dim="year")
    array = array.ffill(dim="year")

    return array


def get_vehicle_fleet_composition(model, vehicle_type) -> Union[xr.DataArray, None]:
    """
    Read the fleet composition csv file and return an `xarray` with dimensions:
    "region", "year", "powertrain", "construction_year", "size"
    :param model: the model to get the fleet composition for
    :param vehicle_type: the type of vehicle to get the fleet composition for
    :return: a multidimensional array with fleet composition data
    """

    if not FILEPATH_FLEET_COMP.is_file():
        raise FileNotFoundError("The fleet composition file could not be found.")

    if model == "remind":
        dataframe = pd.read_csv(
            FILEPATH_FLEET_COMP, sep=get_delimiter(filepath=FILEPATH_FLEET_COMP)
        )
    elif model == "image":
        dataframe = pd.read_csv(
            FILEPATH_IMAGE_TRUCKS_FLEET_COMP,
            sep=get_delimiter(filepath=FILEPATH_FLEET_COMP),
        )
    else:
        return None

    dataframe = dataframe.loc[~dataframe["region"].isnull()]

    with open(VEHICLES_MAP, "r", encoding="utf-8") as stream:
        size_ftr = yaml.safe_load(stream)[vehicle_type]["sizes"]

    dataframe = dataframe.loc[dataframe["size"].isin(size_ftr)]

    if len(dataframe) > 0:
        arr = (
            dataframe.groupby(
                ["region", "year", "powertrain", "construction_year", "size"]
            )
            .sum()["vintage_demand_vkm"]
            .to_xarray()
        )
        arr = arr.fillna(0)

        return arr

    return None


def fix_efficiencies(data: xr.DataArray, min_year: int) -> xr.DataArray:
    """
    Fix the efficiency data to ensure plausibility.

    """

    # If we are looking at a year post 2020
    # and the ratio in efficiency change is inferior to 1
    # we correct it to 1, as we do not accept
    # that efficiency degrades over time
    data.loc[dict(year=[y for y in data.year.values if y > 2020])] = np.clip(
        data.loc[dict(year=[y for y in data.year.values if y > 2020])],
        1,
        None,
    )

    # Inversely, if we are looking at a year prior to 2020
    # and the ratio in efficiency change is superior to 1
    # we correct it to 1, as we do not accept
    # that efficiency in the past was higher than now
    data.loc[dict(year=[y for y in data.year.values if y < 2020])] = np.clip(
        data.loc[dict(year=[y for y in data.year.values if y < 2020])],
        None,
        1,
    )

    # ensure that efficiency can not decrease over time
    while data.diff(dim="year").min().values < 0:
        diff = data.diff(dim="year")
        diff = xr.concat([data.sel(year=min_year), diff], dim="year")
        diff.values[diff.values > 0] = 0
        diff *= -1
        data += diff

    # convert NaNs
    # back-fill missing values with nearest available
    data = data.fillna(data.ffill(dim="variables"))
    # forward-fill missing values with nearest available
    data = data.fillna(data.bfill(dim="variables"))
    data = data.fillna(1)

    # only consider efficiency change between
    # 50% and 300% relative to 2020
    data.values = np.clip(data, 0.5, None)
    data.values = np.clip(data, None, 3)

    return data


def flatten(list_to_flatten):
    rt = []
    for i in list_to_flatten:
        if isinstance(i, list):
            rt.extend(flatten(i))
        else:
            rt.append(i)
    return rt


class IAMDataCollection:
    """
    :var model: name of the IAM model (e.g., "remind")
    :var pathway: name of the IAM scenario (e.g., "SSP2-Base")
    :var year: year to produce the database for
    :var system_model: "cutoff" or "consequential".
    """

    def __init__(
        self,
        model: str,
        pathway: str,
        year: int,
        filepath_iam_files: Path,
        key: bytes,
        external_scenarios: dict = None,
        system_model: str = "cutoff",
        system_model_args: dict = None,
        gains_scenario: str = "CLE",
        use_absolute_efficiency: bool = False,
    ) -> None:
        self.model = model
        self.pathway = pathway
        self.year = year
        self.external_scenarios = external_scenarios
        self.system_model_args = system_model_args
        self.use_absolute_efficiency = use_absolute_efficiency
        self.min_year = 2005
        self.max_year = 2100
        key = key or None

        electricity_prod_vars = self.__get_iam_variable_labels(
            IAM_ELEC_VARS, variable="iam_aliases"
        )
        electricity_eff_vars = self.__get_iam_variable_labels(
            IAM_ELEC_VARS, variable="eff_aliases"
        )

        fuel_prod_vars = self.__get_iam_variable_labels(
            IAM_FUELS_VARS, variable="iam_aliases"
        )
        fuel_eff_vars = self.__get_iam_variable_labels(
            IAM_FUELS_VARS, variable="eff_aliases"
        )

        cement_prod_vars = self.__get_iam_variable_labels(
            IAM_CEMENT_VARS, variable="iam_aliases"
        )

        cement_energy_vars = self.__get_iam_variable_labels(
            IAM_CEMENT_VARS, variable="energy_use_aliases"
        )

        cement_eff_vars = self.__get_iam_variable_labels(
            IAM_CEMENT_VARS, variable="eff_aliases"
        )

        steel_prod_vars = self.__get_iam_variable_labels(
            IAM_STEEL_VARS, variable="iam_aliases"
        )

        steel_energy_vars = self.__get_iam_variable_labels(
            IAM_STEEL_VARS, variable="energy_use_aliases"
        )

        steel_eff_vars = self.__get_iam_variable_labels(
            IAM_STEEL_VARS, variable="eff_aliases"
        )

        dac_prod_vars = self.__get_iam_variable_labels(
            IAM_DAC_VARS, variable="iam_aliases"
        )

        dac_heat_vars = self.__get_iam_variable_labels(
            IAM_DAC_VARS, variable="heat_use_aliases"
        )
        dac_electricity_vars = self.__get_iam_variable_labels(
            IAM_DAC_VARS, variable="electricity_use_aliases"
        )

        biomass_prod_vars = self.__get_iam_variable_labels(
            IAM_BIOMASS_VARS, variable="iam_aliases"
        )

        biomass_eff_vars = self.__get_iam_variable_labels(
            IAM_BIOMASS_VARS, variable="eff_aliases"
        )

        land_use_vars = self.__get_iam_variable_labels(
            IAM_CROPS_VARS, variable="land_use"
        )
        land_use_change_vars = self.__get_iam_variable_labels(
            IAM_CROPS_VARS, variable="land_use_change"
        )

        carbon_capture_vars = self.__get_iam_variable_labels(
            IAM_CARBON_CAPTURE_VARS, variable="iam_aliases"
        )

        other_vars = self.__get_iam_variable_labels(
            IAM_OTHER_VARS, variable="iam_aliases"
        )

        # new_vars is a list of all variables that are declared above

        new_vars = (
            list(electricity_prod_vars.values())
            + list(electricity_eff_vars.values())
            + list(fuel_prod_vars.values())
            + list(fuel_eff_vars.values())
            + list(cement_prod_vars.values())
            + list(cement_energy_vars.values())
            + list(cement_eff_vars.values())
            + list(steel_prod_vars.values())
            + list(steel_energy_vars.values())
            + list(dac_prod_vars.values())
            + list(biomass_prod_vars.values())
            + list(biomass_eff_vars.values())
            + list(land_use_vars.values())
            + list(land_use_change_vars.values())
            + list(carbon_capture_vars.values())
            + list(other_vars.values())
        )

        # flatten the list of lists
        new_vars = flatten(new_vars)

        # if "liquid fossil fuels" is in the list of fuel variables
        # we add the split of gasoline, diesel, LPG and kerosene
        # to `data`, because it means it's not already in the IAM file.
        data = self.__get_iam_data(
            key=key,
            filedir=filepath_iam_files,
            variables=new_vars,
            split_fossil_liquid_fuels=(
                fuel_prod_vars if "liquid fossil fuels" in fuel_prod_vars else None
            ),
        )

        self.data = data

        self.regions = data.region.values.tolist()
        self.system_model = system_model

        self.gains_data_EU = get_gains_EU_data()
        self.gains_data_IAM = get_gains_IAM_data(
            self.model, gains_scenario=gains_scenario
        )

        self.electricity_markets = self.__fetch_market_data(
            data=data,
            input_vars=electricity_prod_vars,
            system_model=self.system_model,
            sector="electricity",
        )

        self.petrol_markets = self.__fetch_market_data(
            data=data,
            input_vars={
                k: v
                for k, v in fuel_prod_vars.items()
                if any(
                    k.lower().startswith(x)
                    for x in [
                        "gasoline",
                        "ethanol",
                        "methanol",
                        "bioethanol",
                        "petrol,",
                    ]
                )
            },
            system_model=self.system_model,
            sector="petrol",
        )

        self.diesel_markets = self.__fetch_market_data(
            data=data,
            input_vars={
                k: v
                for k, v in fuel_prod_vars.items()
                if any(
                    k.lower().startswith(x)
                    for x in [
                        "diesel",
                        "biodiesel",
                    ]
                )
            },
            system_model=self.system_model,
            sector="diesel",
        )

        self.gas_markets = self.__fetch_market_data(
            data=data,
            input_vars={
                k: v
                for k, v in fuel_prod_vars.items()
                if any(
                    k.lower().startswith(x)
                    for x in ["biogas", "methane", "natural gas", "biomethane"]
                )
            },
            system_model=self.system_model,
            sector="gas",
        )

        self.hydrogen_markets = self.__fetch_market_data(
            data=data,
            input_vars={
                k: v
                for k, v in fuel_prod_vars.items()
                if any(
                    k.lower().startswith(x)
                    for x in [
                        "hydrogen",
                    ]
                )
            },
            system_model=self.system_model,
            sector="hydrogen",
        )

        self.kerosene_markets = self.__fetch_market_data(
            data=data,
            input_vars={
                k: v
                for k, v in fuel_prod_vars.items()
                if any(
                    k.lower().startswith(x)
                    for x in [
                        "kerosene",
                    ]
                )
            },
            system_model=self.system_model,
            sector="kerosene",
        )

        self.lpg_markets = self.__fetch_market_data(
            data=data,
            input_vars={
                k: v
                for k, v in fuel_prod_vars.items()
                if any(
                    k.lower().startswith(x)
                    for x in [
                        "liquefied petroleum gas",
                    ]
                )
            },
            system_model=self.system_model,
            sector="lpg",
        )

        self.cement_markets = self.__fetch_market_data(
            data=data,
            input_vars=cement_prod_vars,
            system_model="cutoff",
            sector="cement",
        )
        self.steel_markets = self.__fetch_market_data(
            data=data, input_vars=steel_prod_vars, system_model="cutoff", sector="steel"
        )
        self.dac_markets = self.__fetch_market_data(
            data=data, input_vars=dac_prod_vars, system_model="cutoff", sector="dac"
        )
        self.biomass_markets = self.__fetch_market_data(
            data=data,
            input_vars=biomass_prod_vars,
            system_model="cutoff",
            sector="biomass",
        )

        self.carbon_capture_rate = self.__get_carbon_capture_rate(
            dict_vars=self.__get_iam_variable_labels(
                IAM_CARBON_CAPTURE_VARS,
                variable="iam_aliases",
            ),
            data=data,
        )

        self.other_vars = self.__fetch_market_data(
            data=data,
            input_vars=other_vars,
            normalize=False,
            system_model="cutoff",
        )

        self.electricity_efficiencies = self.get_iam_efficiencies(
            data=data,
            efficiency_labels=electricity_eff_vars,
            use_absolute_efficiency=self.use_absolute_efficiency,
        )

        self.cement_efficiencies = self.get_iam_efficiencies(
            data=data,
            efficiency_labels=cement_eff_vars,
            energy_labels=cement_energy_vars,
            production_labels=cement_prod_vars,
        )

        self.steel_efficiencies = self.get_iam_efficiencies(
            data=data,
            production_labels=steel_prod_vars,
            energy_labels=steel_energy_vars,
            efficiency_labels=steel_eff_vars,
        )

        self.petrol_efficiencies = self.get_iam_efficiencies(
            data=data,
            efficiency_labels={
                k: v
                for k, v in fuel_eff_vars.items()
                if any(
                    k.lower().startswith(x)
                    for x in ["gasoline", "ethanol", "methanol", "bioethanol"]
                )
            },
        )

        self.diesel_efficiencies = self.get_iam_efficiencies(
            data=data,
            efficiency_labels={
                k: v
                for k, v in fuel_eff_vars.items()
                if any(
                    k.lower().startswith(x)
                    for x in [
                        "diesel",
                        "biodiesel",
                    ]
                )
            },
        )

        self.gas_efficiencies = self.get_iam_efficiencies(
            data=data,
            efficiency_labels={
                k: v
                for k, v in fuel_eff_vars.items()
                if any(
                    k.lower().startswith(x)
                    for x in ["biogas", "methane", "natural gas", "biomethane"]
                )
            },
        )

        self.hydrogen_efficiencies = self.get_iam_efficiencies(
            data=data,
            efficiency_labels={
                k: v
                for k, v in fuel_eff_vars.items()
                if any(
                    k.lower().startswith(x)
                    for x in [
                        "hydrogen",
                    ]
                )
            },
        )

        self.kerosene_efficiencies = self.get_iam_efficiencies(
            data=data,
            efficiency_labels={
                k: v
                for k, v in fuel_eff_vars.items()
                if any(
                    k.lower().startswith(x)
                    for x in [
                        "kerosene",
                    ]
                )
            },
        )

        self.lpg_efficiencies = self.get_iam_efficiencies(
            data=data,
            efficiency_labels={
                k: v
                for k, v in fuel_eff_vars.items()
                if any(
                    k.lower().startswith(x)
                    for x in [
                        "liquefied petroleum gas",
                    ]
                )
            },
        )

        self.dac_heat_efficiencies = self.get_iam_efficiencies(
            data=data,
            production_labels=dac_prod_vars,
            energy_labels=dac_heat_vars,
        )

        self.dac_electricity_efficiencies = self.get_iam_efficiencies(
            data=data,
            production_labels=dac_prod_vars,
            energy_labels=dac_electricity_vars,
        )

        self.land_use = self.__get_iam_production_volumes(
            data=data, input_vars=land_use_vars, fill=True
        )

        self.land_use_change = self.__get_iam_production_volumes(
            data=data, input_vars=land_use_change_vars, fill=True
        )

        self.trsp_cars = get_vehicle_fleet_composition(self.model, vehicle_type="car")

        self.trsp_trucks = get_vehicle_fleet_composition(
            self.model, vehicle_type="truck"
        )

        self.trsp_buses = get_vehicle_fleet_composition(self.model, vehicle_type="bus")

        self.metals = get_metals_intensity_factors_data()

        self.production_volumes = self.__get_iam_production_volumes(
            data=data,
            input_vars={
                **electricity_prod_vars,
                **fuel_prod_vars,
                **cement_prod_vars,
                **steel_prod_vars,
                **dac_prod_vars,
                **biomass_prod_vars,
            },
        )

        self.coal_power_plants = self.fetch_external_data_coal_power_plants()

    def __get_iam_variable_labels(
        self, filepath: Path, variable: str
    ) -> Dict[str, Union[str, List[str]]]:
        """
        Loads a csv file into a dictionary.
        This dictionary contains common terminology to ``premise``
        (fuel names, electricity production technologies, etc.) and its
        equivalent variable name in the IAM file.
        :return: dictionary that contains fuel production names equivalence
        """

        dict_vars = {}

        with open(filepath, "r", encoding="utf-8") as stream:
            out = yaml.safe_load(stream)

        for key, values in out.items():
            if variable in values:
                if self.model in values[variable]:
                    if values[variable][self.model] is not None:
                        dict_vars[key] = values[variable][self.model]

        return dict_vars

    def __get_iam_data(
        self,
        key: bytes,
        filedir: Path,
        variables: List,
        split_fossil_liquid_fuels: dict = None,
    ) -> xr.DataArray:
        """
        Read the IAM result file and return an `xarray` with dimensions:

        * region
        * variable
        * year

        :param key: encryption key, if provided by user
        :param filedir: file path to IAM file
        :param variables: list of variables to extract from IAM file

        :return: a multidimensional array with IAM data

        """

        # find file in directory which name contains both self.model and self.pathway
        # Walk through the directory
        filepath = ""
        for root, dirs, files in os.walk(filedir):
            for file in files:
                # Check if both model and pathway are present in the filename
                if self.model in file and self.pathway in file:
                    filepath = Path(os.path.join(root, file))

        if filepath == "":
            raise FileNotFoundError(
                f"Could not find any file containing both {self.model} and {self.pathway} in {filedir}"
            )

        if key is None:
            # Uses a non-encrypted file
            # if extension is ".csv"
            if filepath.suffix in [".csv", ".mif"]:
                print(f"Reading {filepath} as csv file")
                with open(filepath, "rb") as file:
                    # read the encrypted data
                    encrypted_data = file.read()
                    # create a temp csv-like file to pass to pandas.read_csv()
                    data = StringIO(str(encrypted_data, "latin-1"))

            elif filepath.suffix in [".xls", ".xlsx"]:
                print(f"Reading {filepath} as excel file")
                data = pd.read_excel(filepath)

            else:
                raise ValueError(
                    f"Extension {filepath.suffix} is not supported. Please use .csv, .mif, .xls or .xlsx."
                )
        else:
            # Uses an encrypted file
            fernet_obj = Fernet(key)
            with open(filepath, "rb") as file:
                # read the encrypted data
                encrypted_data = file.read()

            # decrypt data
            decrypted_data = fernet_obj.decrypt(encrypted_data)
            data = StringIO(str(decrypted_data, "latin-1"))

        if filepath.suffix in [".csv", ".mif"]:
            dataframe = pd.read_csv(
                data,
                sep=get_delimiter(data=copy.copy(data).readline()),
                encoding="latin-1",
            )
        else:
            dataframe = data

        # if a column name can be an integer
        # we convert it to an integer
        new_cols = {c: int(c) if str(c).isdigit() else c for c in dataframe.columns}
        dataframe = dataframe.rename(columns=new_cols)

        # remove any column that is a string
        # and that is not any of "Region", "Variable", "Unit"
        for col in dataframe.columns:
            if isinstance(col, str):
                if col.lower() not in ["region", "variable", "unit"]:
                    dataframe = dataframe.drop(col, axis=1)

        # identify the lowest and highest column name that is numeric
        # and consider it the minimum year
        self.min_year = min(x for x in dataframe.columns if isinstance(x, int))
        self.max_year = max(x for x in dataframe.columns if isinstance(x, int))

        dataframe = dataframe.reset_index()

        # remove "index" column
        if "index" in dataframe.columns:
            dataframe = dataframe.drop("index", axis=1)

        # convert all column names that are string to lower case
        dataframe.columns = [
            x.lower() if isinstance(x, str) else x for x in dataframe.columns
        ]

        # if split_fossil_liquid_fuels is not None
        # we add the split of gasoline, diesel, LPG and kerosene

        if split_fossil_liquid_fuels is not None:
            # get the split of gasoline, diesel, LPG and kerosene
            df = get_oil_product_volumes(self.model)
            variable_liquid_fuel = split_fossil_liquid_fuels["liquid fossil fuels"]

            for fuel_var, iam_var in split_fossil_liquid_fuels.items():
                if iam_var not in dataframe["variable"].unique():
                    new_fuel_df = dataframe.loc[
                        dataframe["variable"] == variable_liquid_fuel
                    ]
                    fuel_share = df[fuel_var].reindex(new_fuel_df["region"])
                    new_fuel_df.loc[:, "variable"] = iam_var

                    cols = [c for c in new_fuel_df.columns if isinstance(c, int)]
                    with pd.option_context("mode.copy_on_write", True):
                        new_fuel_df.loc[:, cols] = (
                            new_fuel_df.loc[:, cols].astype(float)
                            * fuel_share.values[:, np.newaxis]
                        )
                    dataframe = pd.concat([dataframe, new_fuel_df])

        # filter out unused variables
        dataframe = dataframe.loc[dataframe["variable"].isin(variables)]

        dataframe = dataframe.rename(columns={"variable": "variables"})

        # make a list of headers that are integer
        headers = [x for x in dataframe.columns if isinstance(x, int)]

        # convert the values in these columns to numeric
        dataframe[headers] = dataframe[headers].apply(pd.to_numeric, errors="coerce")

        array = (
            dataframe.melt(
                id_vars=["region", "variables", "unit"],
                var_name="year",
                value_name="value",
            )[["region", "variables", "year", "unit", "value"]]
            .groupby(["region", "variables", "year"])["value"]
            .mean()
            .to_xarray()
        )

        # add the unit as an attribute, as a dictionary with variables as keys
        array.attrs["unit"] = dict(
            dataframe.groupby("variables")["unit"].first().to_dict().items()
        )

        return array

    def __fetch_market_data(
        self,
        data: xr.DataArray,
        input_vars: dict,
        system_model: str,
        normalize: bool = True,
        sector: str = None,
    ) -> [xr.DataArray, None]:
        """
        This method retrieves the market share for each technology,
        for a specified year, for each region provided by the IAM.

        :return: a multidimensional array with technologies market share
        for a given year, for all regions.
        """

        # Check if the year specified is within the range of years given by the IAM
        assert (
            data.year.values.min() <= self.year <= data.year.values.max()
        ), f"{self.year} is outside of the boundaries of the IAM file: {data.year.values.min()}-{data.year.values.max()}"

        # check if values of input_vars are strings or lists
        if any(isinstance(x, list) for x in input_vars.values()):
            vars = list(chain.from_iterable(input_vars.values()))
        else:
            vars = list(input_vars.values())

        missing_vars = set(vars) - set(data.variables.values)

        if missing_vars:
            print_missing_variables(missing_vars)

        available_vars = list(set(vars) - missing_vars)

        if available_vars:
            market_data = data.loc[:, available_vars, :]
        else:
            return None

        if any(isinstance(x, list) for x in input_vars.values()):
            rev_input_vars = {x: k for k, v in input_vars.items() for x in v}
        else:
            rev_input_vars = {v: k for k, v in input_vars.items()}

        market_data.coords["variables"] = [
            rev_input_vars[v] for v in market_data.variables.values
        ]

        # if duplicates in market_data.coords["variables"]
        # we sum them
        if len(market_data.coords["variables"].values.tolist()) != len(
            set(market_data.coords["variables"].values.tolist())
        ):
            market_data = market_data.groupby("variables").sum(dim="variables")

        if system_model == "consequential":
            market_data = consequential_method(
                market_data, self.year, self.system_model_args, sector
            )
        else:
            if normalize is True:
                market_data /= market_data.groupby("region").sum(dim="variables")

        # back-fill nans
        market_data = market_data.bfill(dim="year")
        # fill NaNs with zeros
        market_data = market_data.fillna(0)

        return market_data

    def get_iam_efficiencies(
        self,
        data: xr.DataArray,
        efficiency_labels: dict = None,
        production_labels: dict = None,
        energy_labels: dict = None,
        use_absolute_efficiency: bool = False,
    ) -> [xr.DataArray, None]:
        """
        This method retrieves efficiency values for the specified sector,
        for a specified year, for each region provided by the IAM.

        :param data: The data to process.
        :param efficiency_labels: The efficiency labels to use.
        :param production_labels: The production labels to use.
        :param energy_labels: The energy labels to use.

        :return: a multidimensional array with sector's technologies market
        share for a given year, for all regions.

        """

        efficiency_labels = efficiency_labels or {}
        production_labels = production_labels or {}
        energy_labels = energy_labels or {}

        # Check if the year specified is within the range of years given by the IAM
        # If the year specified is not contained within the range of years given by the IAM
        if self.year < data.year.values.min() or self.year > data.year.values.max():
            raise KeyError(
                f"{self.year} is outside of the boundaries "
                f"of the IAM file: {data.year.values.min()}-{data.year.values.max()}"
            )

        # Finally, if the specified year falls in between two periods provided by the IAM
        # Interpolation between two periods
        if efficiency_labels:
            missing_vars = set(efficiency_labels.values()) - set(data.variables.values)
            if missing_vars:
                print_missing_variables(missing_vars)

            available_vars = list(set(efficiency_labels.values()) - missing_vars)
            rev_eff_labels = {v: k for k, v in efficiency_labels.items()}

            if available_vars:
                eff_data = data.loc[dict(variables=available_vars)]
                eff_data.coords["variables"] = [
                    rev_eff_labels[x] for x in eff_data.variables.values
                ]
                # convert zero values to nan
                # and back-fill and forward-fill missing values with nearest available
                eff_data = eff_data.where(eff_data != 0)
                eff_data = eff_data.bfill(dim="year")
                eff_data = eff_data.ffill(dim="year")

            else:
                return None

        elif production_labels and energy_labels:
            eff_data = xr.DataArray(dims=["variables"], coords={"variables": []})
            for k, v in production_labels.items():

                # check that each element of energy.values() is in data.variables.values
                # knowing that energy.values() is a list of lists
                # and that each element of prod.values() is in data.variables.values
                _ = lambda x: (
                    x
                    if isinstance(x, list)
                    else [
                        x,
                    ]
                )

                if all(
                    var in data.variables.values for var in energy_labels[k]
                ) and all(x in data.variables.values for x in _(v)):
                    if isinstance(v, list):
                        d = data.loc[:, energy_labels[k], :].sum(
                            dim="variables"
                        ) / data.loc[:, list(v), :].sum(dim="variables")
                        # add dimension "variables" to d
                        d = d.expand_dims(dim="variables")
                        # add a coordinate "variables" to d
                        d.coords["variables"] = [
                            k,
                        ]
                    else:
                        d = (
                            data.loc[:, energy_labels[k], :].sum(dim="variables")
                            / data.loc[:, v, :]
                        )
                    # convert inf to Nan
                    d = d.where(d != np.inf)
                    # back-fill nans
                    d = d.bfill(dim="year")
                    # forward-fill nans
                    d = d.ffill(dim="year")

                else:
                    # fill d with ones
                    d = xr.ones_like(data.loc[:, data.variables[0], :])

                eff_data = xr.concat([eff_data, d], dim="variables")

            eff_data.coords["variables"] = list(production_labels.keys())
        else:
            return None

        if use_absolute_efficiency is False:
            # efficiency expressed
            eff_data /= eff_data.sel(year=2020)

            if len(efficiency_labels) == 0 or any(
                "specific" in x.lower() for x in efficiency_labels.values()
            ):
                # we are dealing with specific energy consumption, not efficiencies
                # we need to convert them to efficiencies
                eff_data = 1 / eff_data

            # fix efficiencies
            eff_data = fix_efficiencies(eff_data, self.min_year)

        else:
            # if absolute efficiencies are used, we need to make sure that
            # the efficiency is not greater than 1
            # otherwise it means they are given as percentages

            # check if any efficiency is greater than 1
            if (eff_data > 1).any():
                # if yes, divide by 100
                eff_data /= 100

            # now check that all efficiencies are between 0 and 1
            eff_data = xr.where(eff_data > 1, 1, eff_data)
            eff_data = xr.where(eff_data < 0, 0, eff_data)

        return eff_data

    def __get_carbon_capture_rate(
        self, dict_vars: Dict[str, str], data: xr.DataArray
    ) -> xr.DataArray:
        """
        Returns a xarray with carbon capture rates for steel and cement production.

        :param dict_vars: dictionary that contains AIM variables to search for
        :param data: IAM data
        :return: a xarray with carbon capture rates, for each year and region
        """

        # If the year specified is not contained within the range of years given by the IAM
        if self.year < data.year.values.min() or self.year > data.year.values.max():
            raise KeyError(
                f"{self.year} is outside of the boundaries "
                f"of the IAM file: {data.year.values.min()}-{data.year.values.max()}"
            )

        # Finally, if the specified year falls in between two periods provided by the IAM
        # Interpolation between two periods

        # if variable is missing, we assume that the rate is 0
        # and that none of the  CO2 emissions are captured

        if isinstance(dict_vars.get("cement - cco2", []), str):
            dict_vars["cement - cco2"] = [
                dict_vars["cement - cco2"],
            ]

        if not any(
            x in data.variables.values.tolist()
            for x in dict_vars.get("cement - cco2", [])
        ):
            print("Cannot find variables for cement capture rate.")
            cement_rate = xr.DataArray(
                np.zeros((len(data.region), len(data.year))),
                coords=[data.region, data.year],
                dims=["region", "year"],
            )

        else:
            cement_rate = data.loc[:, dict_vars["cement - cco2"], :].sum(
                dim=["variables"]
            ) / data.loc[:, dict_vars["cement - co2"], :].sum(dim=["variables"])

        cement_rate.coords["variables"] = "cement"

        if isinstance(dict_vars.get("steel - cco2", []), str):
            dict_vars["steel - cco2"] = [
                dict_vars["steel - cco2"],
            ]

        if not any(
            x in data.variables.values.tolist()
            for x in dict_vars.get("steel - cco2", [])
        ):
            print("Cannot find variables for steel capture rate.")
            steel_rate = xr.DataArray(
                np.zeros((len(data.region), len(data.year))),
                coords=[data.region, data.year],
                dims=["region", "year"],
            )
        else:
            steel_rate = data.loc[:, dict_vars["steel - cco2"], :].sum(
                dim="variables"
            ) / data.loc[:, dict_vars["steel - co2"], :].sum(dim="variables")

        steel_rate.coords["variables"] = "steel"

        rate = xr.concat([cement_rate, steel_rate], dim="variables")

        # forward fill missing values
        rate = rate.ffill(dim="year")

        rate = rate.fillna(0)

        # we need to fix the rate for "World"
        # as it is sometimes neglected in the
        # IAM files

        if "World" in rate.region.values.tolist():
            if not any(
                x in data.variables.values.tolist()
                for x in dict_vars.get("cement - cco2", [])
            ):
                rate.loc[dict(region="World", variables="cement")] = 0
            else:
                try:
                    rate.loc[dict(region="World", variables="cement")] = (
                        data.loc[
                            dict(
                                region=[r for r in self.regions if r != "World"],
                                variables=dict_vars["cement - cco2"],
                            )
                        ]
                        .sum(dim=["variables", "region"])
                        .values
                        / data.loc[
                            dict(
                                region=[r for r in self.regions if r != "World"],
                                variables=dict_vars["cement - co2"],
                            )
                        ]
                        .sum(dim=["variables", "region"])
                        .values
                    )

                except ZeroDivisionError:
                    rate.loc[dict(region="World", variables="cement")] = 0

                try:
                    rate.loc[dict(region="World", variables="steel")] = data.loc[
                        dict(
                            region=[r for r in self.regions if r != "World"],
                            variables=dict_vars["steel - cco2"],
                        )
                    ].sum(dim=["variables", "region"]) / data.loc[
                        dict(
                            region=[r for r in self.regions if r != "World"],
                            variables=dict_vars["steel - co2"],
                        )
                    ].sum(
                        dim=["variables", "region"]
                    )

                except ZeroDivisionError:
                    rate.loc[dict(region="World", variables="steel")] = 0

            if not any(
                x in data.variables.values.tolist()
                for x in dict_vars.get("steel - cco2", [])
            ):
                rate.loc[dict(region="World", variables="steel")] = 0
            else:
                rate.loc[dict(region="World", variables="steel")] = (
                    data.loc[
                        dict(
                            region=[r for r in self.regions if r != "World"],
                            variables=dict_vars["steel - cco2"],
                        )
                    ]
                    .sum(dim=["variables", "region"])
                    .values
                    / data.loc[
                        dict(
                            region=[r for r in self.regions if r != "World"],
                            variables=dict_vars["steel - co2"],
                        )
                    ]
                    .sum(dim=["variables", "region"])
                    .values
                )

        # we ensure that the rate can only be between 0 and 1
        rate.values = np.clip(rate, 0, 1)

        # values under 0.001 are considered as 0
        rate = xr.where(rate < 0.001, 0, rate)

        return rate

    def __get_iam_production_volumes(
        self, input_vars, data, fill: bool = False
    ) -> [xr.DataArray, None]:
        """
        Returns n xarray with production volumes for different sectors:
        electricity, steel, cement, fuels.
        This is used to build markets: we use
        the production volumes of each region for example,
        to build the World market.
        :param input_vars: a dictionary that contains
        common labels as keys, and IAM labels as values.
        :param data: IAM data
        :return: a xarray with production volumes for
        different commodities (electricity, cement, etc.)
        """

        def flatten_list_to_strings(input_list):
            result = []
            for element in input_list:
                if isinstance(element, str):
                    result.append(element)
                elif isinstance(element, list):
                    result.extend(flatten_list_to_strings(element))
            return result

        # If the year specified is not contained within the range of years given by the IAM
        if self.year < data.year.values.min() or self.year > data.year.values.max():
            raise KeyError(
                f"{self.year} is outside of the boundaries "
                f"of the IAM file: {data.year.values.min()}-{data.year.values.max()}"
            )

        vars = flatten_list_to_strings(input_vars.values())

        missing_vars = set(vars) - set(data.variables.values)

        if missing_vars:
            print_missing_variables(missing_vars)

        available_vars = list(set(vars) - missing_vars)

        if available_vars:
            data_to_return = data.loc[:, available_vars, :]
        else:
            return None

        if any(isinstance(x, list) for x in input_vars.values()):
            rev_input_vars = {}
            for k, v in input_vars.items():
                if isinstance(v, list):
                    for x in v:
                        rev_input_vars[x] = k
                else:
                    rev_input_vars[v] = k
        else:
            rev_input_vars = {v: k for k, v in input_vars.items()}

        data_to_return.coords["variables"] = [
            rev_input_vars[v] for v in data_to_return.variables.values
        ]

        # if duplicates in market_data.coords["variables"]
        # we sum them
        if len(data_to_return.coords["variables"].values.tolist()) != len(
            set(data_to_return.coords["variables"].values.tolist())
        ):
            data_to_return = data_to_return.groupby("variables").sum(dim="variables")

        if fill:
            # if fill, we fill zero values
            # with the nearest year's value
            # first, convert zero values to NaNs
            data_to_return = data_to_return.where(data_to_return != 0)
            # then, backfill
            data_to_return = data_to_return.bfill(dim="year")
            # then, forward fill
            data_to_return = data_to_return.ffill(dim="year")

        return data_to_return

    def get_external_data(self, external_scenarios: list):
        """
        Fetch data from external sources.
        :param external_scenarios: a list of dictionaries
        with keys "scenario" and "data"
        :return: a dictionary with data

        """
        data = {}

        for i, external_scenario in enumerate(external_scenarios):
            scenario, dp = external_scenario["scenario"], external_scenario["data"]
            data[i] = {}

            resource = dp.get_resource("scenario_data")
            # getting scenario data in binary format
            scenario_data = resource.raw_read()
            df = pd.read_csv(
                BytesIO(scenario_data),
            )
            # set headers from first row
            try:
                df.columns = resource.headers
            except ValueError as err:
                raise ValueError(
                    f"The number of headers in scenario data file are not correct. {err}"
                    f"Check that the values in the scenario data file are separated by commas, not semicolons."
                ) from err

            resource = dp.get_resource("config")
            config_file = yaml.safe_load(resource.raw_read())

            if "production pathways" in config_file:
                variables = {}
                for k, v in config_file["production pathways"].items():
                    try:
                        variables[k] = v["production volume"]["variable"]
                    except KeyError:
                        continue

                subset = df.loc[
                    (df["scenario"] == scenario)
                    & (df["variables"].isin(variables.values())),
                    "region":,
                ]

                array = (
                    subset.melt(
                        id_vars=["region", "variables", "unit"],
                        var_name="year",
                        value_name="value",
                    )[["region", "variables", "year", "value"]]
                    .groupby(["region", "variables", "year"])["value"]
                    .mean()
                    .to_xarray()
                )

                # convert to float64
                array = array.astype(np.float64)
                # convert year dim to int64
                array.coords["year"] = array.coords["year"].astype(np.int64)

                data[i]["production volume"] = array
                regions = subset["region"].unique().tolist()
                data[i]["regions"] = regions

                variables = {}
                if "production pathways" in config_file:
                    for k, v in config_file["production pathways"].items():
                        try:
                            variables[k] = [e["variable"] for e in v["efficiency"]]
                        except KeyError:
                            continue

                if "markets" in config_file:
                    for m, market in enumerate(config_file["markets"]):
                        try:
                            variables[f"market {m}"] = [
                                e["variable"] for e in market["efficiency"]
                            ]
                        except KeyError:
                            continue

                if len(variables) > 0:
                    subset = df.loc[
                        (df["scenario"] == scenario)
                        & (df["variables"].isin(list(chain(*variables.values())))),
                        "region":,
                    ]

                    array = (
                        subset.melt(
                            id_vars=["region", "variables", "unit"],
                            var_name="year",
                            value_name="value",
                        )[["region", "variables", "year", "value"]]
                        .groupby(["region", "variables", "year"])["value"]
                        .mean()
                        .to_xarray()
                    )
                    # convert to float64
                    array = array.astype(np.float64)
                    # convert year dim to int64
                    array.coords["year"] = array.coords["year"].astype(np.int64)

                    ref_years = {}

                    if "production pathways" in config_file:
                        for v in config_file["production pathways"].values():
                            for e, f in v.items():
                                if e == "efficiency":
                                    for x in f:
                                        ref_years[x["variable"]] = x.get(
                                            "reference year", None
                                        )

                    if "markets" in config_file:
                        for market in config_file["markets"]:
                            for e, f in market.items():
                                if f == "efficiency":
                                    for x in f["efficiency"]:
                                        ref_years[x["variable"]] = x.get(
                                            "reference year", None
                                        )

                    for y, ref_year in ref_years.items():
                        if ref_year is None:
                            # use the earliest year in `array`
                            ref_years[y] = array.year.values.min()

                    for v, y in ref_years.items():
                        array.loc[{"variables": v}] = array.loc[
                            {"variables": v}
                        ] / array.loc[{"variables": v}].sel(year=int(y))

                    # convert NaNs to ones
                    array = array.fillna(1)

                    data[i]["efficiency"] = array

        return data

    def fetch_external_data_coal_power_plants(self):
        """
        Fetch data on coal power plants from external sources.
        Source:
        Oberschelp, C., Pfister, S., Raptis, C.E. et al.
        Global emission hotspots of coal power generation.
        Nat Sustain 2, 113–121 (2019).
        https://doi.org/10.1038/s41893-019-0221-6

        """

        df = pd.read_csv(COAL_POWER_PLANTS_DATA, sep=",", index_col=False)
        # rename columns
        new_cols = {
            "ISO2": "country",
            "NET_ELECTRICITY_GENERATION_MWH": "generation",
            "FUEL_INPUT_LHV_MJ": "fuel input",
            "NET_ELECTRICAL_EFFICIENCY": "efficiency",
            "CHP_PLANT": "CHP",
            "PLANT_FUEL": "fuel",
            "PLANT_EMISSION_CO2_KG": "CO2",
            "PLANT_EMISSION_CH4_KG": "CH4",
            "PLANT_EMISSION_SO2_KG": "SO2",
            "PLANT_EMISSION_NOX_KG": "NOx",
            "PLANT_EMISSION_PM_2.5_KG": "PM <2.5",
            "PLANT_EMISSION_PM_10_TO_2.5_KG": "PM 10 - 2.5",
            "PLANT_EMISSION_PM_GR_10_KG": "PM > 10",
            "PLANT_EMISSION_HG_0_KG": "HG0",
            "PLANT_EMISSION_HG_2P_KG": "HG2",
            "PLANT_EMISSION_HG_P_KG": "HGp",
        }
        df = df.rename(columns=new_cols)

        # drop columns
        df = df.drop(columns=[c for c in df.columns if c not in new_cols.values()])

        # rename Bituminous fuel type as Anthracite
        df.loc[:, "fuel"] = df.loc[:, "fuel"].replace(
            "Bituminous coal", "Anthracite coal"
        )

        # rename Subbituminous  and Coal blend fuel type as Lignite
        df["fuel"] = df["fuel"].replace("Subbituminous coal", "Lignite coal")
        df["fuel"] = df["fuel"].replace("Coal blend", "Lignite coal")

        # convert to xarray
        # with dimensions: country and CHP
        # with variables as avergage: generation, efficiency, CO2, CH4, SO2,
        # NOx, PM <2.5, PM 10 - 2.5, PM > 10, HG0, HG2, HGp
        # and ignore the following variables: fuel input

        df = df.drop(columns=["fuel input"])
        array = (
            df.melt(
                id_vars=[
                    "country",
                    "CHP",
                    "fuel",
                ],
                var_name="variable",
                value_name="value",
            )
            .groupby(["country", "CHP", "fuel", "variable"])["value"]
            .mean()
            .to_xarray()
        )

        return array
