"""
data_collection.py contains the IAMDataCollection class which collects a number of data,
mostly from the IAM file. This class will have offer market shares, efficiency
and emission values for different sectors, carbon capture rates, etc.
Additional external sources of data have to be used as well, notably for cement
production (GNR data), and for non-CO2 emissions (GAINS data).
"""


import csv
from io import StringIO
from itertools import chain
from pathlib import Path
from typing import Dict, List, Union

import numpy as np
import pandas as pd
import xarray as xr
import yaml
from cryptography.fernet import Fernet

from . import DATA_DIR

IAM_ELEC_VARS = DATA_DIR / "electricity" / "electricity_tech_vars.yml"
IAM_FUELS_VARS = DATA_DIR / "fuels" / "fuel_tech_vars.yml"
IAM_BIOMASS_VARS = DATA_DIR / "electricity" / "biomass_vars.yml"
IAM_CROPS_VARS = DATA_DIR / "fuels" / "crops_properties.yml"
IAM_CEMENT_VARS = DATA_DIR / "cement" / "cement_tech_vars.yml"
IAM_STEEL_VARS = DATA_DIR / "steel" / "steel_tech_vars.yml"
IAM_OTHER_VARS = DATA_DIR / "utils" / "report" / "other_vars.yaml"
IAM_LIFETIMES = DATA_DIR / "lifetimes.csv"
FILEPATH_FLEET_COMP = (
    DATA_DIR / "iam_output_files" / "fleet_files" / "fleet_all_vehicles.csv"
)
FILEPATH_IMAGE_TRUCKS_FLEET_COMP = (
    DATA_DIR / "iam_output_files" / "fleet_files" / "image_fleet_trucks.csv"
)
VEHICLES_MAP = DATA_DIR / "transport" / "vehicles_map.yaml"
GAINS_TO_IAM_FILEPATH = DATA_DIR / "GAINS_emission_factors" / "GAINStoREMINDtechmap.csv"
GNR_DATA = DATA_DIR / "cement" / "additional_data_GNR.csv"
IAM_CARBON_CAPTURE_VARS = DATA_DIR / "utils" / "carbon_capture_vars.yml"
CROPS_PROPERTIES = DATA_DIR / "fuels" / "crops_properties.yml"


def get_crops_properties() -> dict:
    """
    Return a dictionary with crop names as keys and IAM labels as values
    relating to land use change CO2 per crop type
    :return: dict
    """
    with open(CROPS_PROPERTIES, "r", encoding="utf-8") as stream:
        crop_props = yaml.safe_load(stream)

    return crop_props


def get_lifetime(list_tech: List) -> np.array:
    """
    Fetch lifetime values for different technologies from a .csv file.
    This is only used for consequential databases.
    :param list_tech: technology labels to find lifetime values for.
    :return: a numpy array with technology lifetime values
    """
    dict_ = {}
    with open(IAM_LIFETIMES, encoding="utf-8") as file:
        reader = csv.reader(file, delimiter=";")
        for row in reader:
            dict_[row[0]] = row[1]

    arr = np.zeros_like(list_tech)

    for i, tech in enumerate(list_tech):
        lifetime = dict_[tech]
        arr[i] = lifetime

    return arr.astype(float)


def get_gnr_data() -> xr.DataArray:
    """
    Read the GNR csv file on cement production and return an `xarray` with dimensions:

    * region
    * year
    * variables
    This data is further used in cement.py.

    :return: a multi-dimensional array with GNR data about cement production

    """
    dataframe = pd.read_csv(GNR_DATA)
    dataframe = dataframe[["region", "year", "variables", "value"]]

    gnr_array = (
        dataframe.groupby(["region", "year", "variables"]).mean()["value"].to_xarray()
    )
    gnr_array = gnr_array.interpolate_na(
        dim="year", method="linear", fill_value="extrapolate"
    )
    gnr_array = gnr_array.interp(year=2020)
    gnr_array = gnr_array.fillna(0)

    return gnr_array


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

    return array


def get_gains_data() -> xr.DataArray:
    """
    Read the GAINS emissions csv file and return an `xarray` with dimensions:

    * region
    * pollutant
    * sector
    * year

    :return: a multi-dimensional array with GAINS emissions data

    """
    filename = "GAINS emission factors.csv"
    filepath = DATA_DIR / "GAINS_emission_factors" / filename

    gains_emi = pd.read_csv(
        filepath,
        skiprows=4,
        names=["year", "region", "GAINS", "pollutant", "pathway", "factor"],
    )
    gains_emi["unit"] = "Mt/TWa"
    gains_emi = gains_emi[gains_emi.pathway == "SSP2"]

    sector_mapping = pd.read_csv(GAINS_TO_IAM_FILEPATH).drop(
        ["noef", "elasticity"], axis=1
    )

    gains_emi = (
        gains_emi.join(sector_mapping.set_index("GAINS"), on="GAINS")
        .dropna()
        .drop(["pathway", "REMIND"], axis=1)
        .pivot_table(
            index=["region", "GAINS", "pollutant", "unit"],
            values="factor",
            columns="year",
        )
    )

    gains_emi = gains_emi.reset_index()
    gains_emi = gains_emi.melt(
        id_vars=["region", "pollutant", "unit", "GAINS"],
        var_name="year",
        value_name="value",
    )[["region", "pollutant", "GAINS", "year", "value"]]
    gains_emi = gains_emi.rename(columns={"GAINS": "sector"})
    array = (
        gains_emi.groupby(["region", "pollutant", "year", "sector"])["value"]
        .mean()
        .to_xarray()
    )

    return array / 8760  # per TWha --> per TWh


def get_vehicle_fleet_composition(model, vehicle_type) -> Union[xr.DataArray, None]:
    """
    Read the fleet composition csv file and return an `xarray` with dimensions:
    "region", "year", "powertrain", "construction_year", "size"
    :param model: the model to get the fleet composition for
    :param vehicle_type: the type of vehicle to get the fleet composition for
    :return: a multi-dimensional array with fleet composition data
    """

    if not FILEPATH_FLEET_COMP.is_file():
        raise FileNotFoundError("The fleet composition file could not be found.")

    if model == "remind":
        dataframe = pd.read_csv(FILEPATH_FLEET_COMP, sep=";")
    else:
        dataframe = pd.read_csv(FILEPATH_IMAGE_TRUCKS_FLEET_COMP, sep=";")

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


class IAMDataCollection:

    """
    :var model: name of the IAM model (e.g., "remind")
    :var pathway: name of the IAM scenario (e.g., "SSP2-Base")
    :var year: year to produce the database for
    :var filepath_iam_files: if a custom file is provided, the filepath to it
    :var key: decryption key, provided by user, if built-in IAM files are to be used
    :var system_model: "attributional" or "consequential" (not yet implemented.
    :var time_horizon: for consequential modelling (not yet implemented.
    """

    def __init__(
        self,
        model: str,
        pathway: str,
        year: int,
        filepath_iam_files: Path,
        key: bytes,
        system_model: str = "attributional",
        time_horizon: int = 30,
        external_scenarios: list = None,
    ) -> None:
        self.model = model
        self.pathway = pathway
        self.year = year
        self.external_scenarios = external_scenarios
        key = key or None

        prod_vars = self.__get_iam_variable_labels(
            IAM_ELEC_VARS, variable="iam_aliases"
        )
        eff_vars = self.__get_iam_variable_labels(IAM_ELEC_VARS, variable="eff_aliases")

        prod_vars.update(
            self.__get_iam_variable_labels(IAM_FUELS_VARS, variable="iam_aliases")
        )
        eff_vars.update(
            self.__get_iam_variable_labels(IAM_FUELS_VARS, variable="eff_aliases")
        )
        prod_vars.update(
            self.__get_iam_variable_labels(IAM_CEMENT_VARS, variable="iam_aliases")
        )
        eff_vars.update(
            self.__get_iam_variable_labels(
                IAM_CEMENT_VARS, variable="energy_use_aliases"
            )
        )
        eff_vars.update(
            self.__get_iam_variable_labels(IAM_CEMENT_VARS, variable="eff_aliases")
        )
        energy_use_vars = self.__get_iam_variable_labels(
            IAM_CEMENT_VARS, variable="energy_use_aliases"
        )
        prod_vars.update(
            self.__get_iam_variable_labels(IAM_STEEL_VARS, variable="iam_aliases")
        )
        eff_vars.update(
            self.__get_iam_variable_labels(
                IAM_STEEL_VARS, variable="energy_use_aliases"
            )
        )
        energy_use_vars.update(
            self.__get_iam_variable_labels(
                IAM_STEEL_VARS, variable="energy_use_aliases"
            )
        )
        prod_vars.update(
            self.__get_iam_variable_labels(IAM_BIOMASS_VARS, variable="iam_aliases")
        )
        eff_vars.update(
            self.__get_iam_variable_labels(IAM_BIOMASS_VARS, variable="eff_aliases")
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

        variables = list(prod_vars.values())
        variables.extend(eff_vars.values())
        variables.extend(energy_use_vars.values())
        variables.extend(land_use_vars.values())
        variables.extend(land_use_change_vars.values())
        variables.extend(carbon_capture_vars.values())
        variables.extend(other_vars.values())
        new_vars = []
        for variable in variables:
            if isinstance(variable, list):
                for sub_var in variable:
                    new_vars.append(sub_var)
            else:
                new_vars.append(variable)

        data = self.__get_iam_data(
            key=key,
            filedir=filepath_iam_files,
            variables=new_vars,
        )

        self.regions = data.region.values.tolist()
        self.system_model = system_model
        self.time_horizon = time_horizon

        gains_data = get_gains_data()
        self.gnr_data = get_gnr_data()

        self.electricity_markets = self.__get_iam_electricity_markets(data=data)
        self.fuel_markets = self.__get_iam_fuel_markets(data=data)

        self.production_volumes = self.__get_iam_production_volumes(
            prod_vars, data=data
        )
        self.carbon_capture_rate = self.__get_carbon_capture_rate(
            dict_vars=self.__get_iam_variable_labels(
                IAM_CARBON_CAPTURE_VARS, variable="iam_aliases"
            ),
            data=data,
        )

        self.other_vars = self.__get_other_iam_vars(data=data)

        electricity_efficiencies = self.__get_iam_electricity_efficiencies(data=data)
        electricity_emissions = self.__get_gains_electricity_emissions(data=gains_data)
        cement_emissions = self.__get_gains_cement_emissions(data=gains_data)
        cement_efficiencies = self.__get_iam_cement_efficiencies(data=data)
        steel_emissions = self.__get_gains_steel_emissions(data=gains_data)
        steel_efficiencies = self.__get_iam_steel_efficiencies(data=data)
        fuel_efficiencies = self.__get_iam_fuel_efficiencies(data=data)

        self.efficiency = xr.concat(
            [
                electricity_efficiencies,
                steel_efficiencies,
                cement_efficiencies,
                fuel_efficiencies,
            ],
            dim="variables",
        )
        self.emissions = xr.concat(
            [
                electricity_emissions,
                steel_emissions,
                cement_emissions,
            ],
            dim="sector",
        )

        if self.model == "image":
            self.land_use = self.__get_iam_land_use(data=data)
            self.land_use_change = self.__get_iam_land_use_change_emissions(data=data)
        else:
            self.land_use = None
            self.land_use_change = None

        self.trsp_cars = get_vehicle_fleet_composition(self.model, vehicle_type="car")
        self.trsp_trucks = get_vehicle_fleet_composition(
            self.model, vehicle_type="truck"
        )
        self.trsp_buses = get_vehicle_fleet_composition(self.model, vehicle_type="bus")

        self.metals = get_metals_intensity_factors_data()

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
                if variable == "gains_aliases":
                    dict_vars[key] = values[variable]
                else:
                    if self.model in values[variable]:
                        dict_vars[key] = values[variable][self.model]

        return dict_vars

    def __get_iam_data(
        self, key: bytes, filedir: Path, variables: List
    ) -> xr.DataArray:
        """
        Read the IAM result file and return an `xarray` with dimensions:

        * region
        * variable
        * year

        :param key: encryption key, if provided by user
        :param filepath: file path to IAM file

        :return: a multi-dimensional array with IAM data

        """

        file_ext = self.model + "_" + self.pathway + ".csv"
        filepath = Path(filedir) / file_ext

        if key is None:
            # Uses a non-encrypted file
            try:
                with open(filepath, "rb") as file:
                    # read the encrypted data
                    encrypted_data = file.read()
            except FileNotFoundError:
                file_ext = self.model + "_" + self.pathway + ".mif"
                filepath = Path(filedir) / file_ext
                with open(filepath, "rb") as file:
                    # read the encrypted data
                    encrypted_data = file.read()

            # create a temp csv-like file to pass to pandas.read_csv()
            data = StringIO(str(encrypted_data, "latin-1"))

        else:
            # Uses an encrypted file
            fernet_obj = Fernet(key)
            with open(filepath, "rb") as file:
                # read the encrypted data
                encrypted_data = file.read()

            # decrypt data
            decrypted_data = fernet_obj.decrypt(encrypted_data)
            data = StringIO(str(decrypted_data, "latin-1"))

        if self.model == "remind":
            dataframe = pd.read_csv(
                data,
                sep=";",
                index_col=["Region", "Variable", "Unit"],
                encoding="latin-1",
            ).drop(columns=["Model", "Scenario"])

            # if new sub-European regions are present, we remove EUR and NEU
            if any(
                x in dataframe.index.get_level_values("Region").unique()
                for x in ["ESC", "DEU", "NEN"]
            ):
                dataframe = dataframe.loc[
                    ~dataframe.index.get_level_values("Region").isin(["EUR", "NEU"])
                ]

            if len(dataframe.columns == 20):
                dataframe.drop(columns=dataframe.columns[-1], inplace=True)

        elif self.model == "image":

            dataframe = pd.read_csv(
                data, index_col=[2, 3, 4], encoding="latin-1", sep=","
            ).drop(columns=["Model", "Scenario"])

        else:
            raise ValueError(
                f"The IAM model name {self.model.upper()} is not valid."
                f"Currently supported: 'REMIND' or 'IMAGE'"
            )

        dataframe.columns = dataframe.columns.astype(int)
        dataframe = dataframe.reset_index()

        dataframe = dataframe.loc[dataframe["Variable"].isin(variables)]

        dataframe = dataframe.rename(
            columns={"Region": "region", "Variable": "variables", "Unit": "unit"}
        )

        array = (
            dataframe.melt(
                id_vars=["region", "variables", "unit"],
                var_name="year",
                value_name="value",
            )[["region", "variables", "year", "value"]]
            .groupby(["region", "variables", "year"])["value"]
            .mean()
            .to_xarray()
        )

        return array

    def __transform_to_marginal_markets(self, data: xr.DataArray) -> xr.DataArray:
        """
        Used for consequential modeling only. Returns marginal market mixes.

        :param data: IAM data
        :return: marginal market mixes
        """

        shape = list(data.shape)
        shape[-1] = 1

        market_shares = xr.DataArray(
            np.zeros(tuple(shape)),
            dims=["region", "variables", "year"],
            coords={
                "region": data.coords["region"],
                "variables": data.variables,
                "year": [self.year],
            },
        )

        for region in data.coords["region"].values:

            current_shares = data.sel(region=region, year=self.year) / data.sel(
                region=region, year=self.year
            ).sum(dim="variables")

            # we first need to calculate the average capital replacement rate of the market
            # which is here defined as the inverse of the production-weighted average lifetime
            lifetime = get_lifetime(current_shares.variables.values.tolist())

            avg_lifetime = np.sum(current_shares.values * lifetime)

            avg_cap_repl_rate = -1 / avg_lifetime

            volume_change = (
                data.sel(region=region)
                .sum(dim="variables")
                .interp(year=self.year + self.time_horizon)
                / data.sel(region=region).sum(dim="variables").interp(year=self.year)
            ) - 1

            # first, we set CHP suppliers to zero
            # as electricity production is not a determining product for CHPs
            tech_to_ignore = ["CHP", "biomethane"]
            data.loc[
                dict(
                    variables=[
                        v
                        for v in data.variables.values
                        if any(x in v for x in tech_to_ignore)
                    ],
                    region=region,
                )
            ] = 0

            # second, we fetch the ratio between production in `self.year` and `self.year` + `time_horizon`
            # for each technology
            market_shares.loc[dict(region=region)] = (
                data.sel(region=region)
                .interp(year=self.year + self.time_horizon)
                .values
                / data.sel(region=region).interp(year=self.year).values
            )[:, None] - 1

            market_shares.loc[dict(region=region)] = market_shares.loc[
                dict(region=region)
            ].round(3)

            # we remove NaNs and np.inf
            market_shares.loc[dict(region=region)].values[
                market_shares.loc[dict(region=region)].values == np.inf
            ] = 0
            market_shares.loc[dict(region=region)] = market_shares.loc[
                dict(region=region)
            ].fillna(0)

            # we fetch the technologies' lifetimes
            lifetime = get_lifetime(market_shares.variables.values)
            # get the capital replacement rate
            # which is here defined as -1 / lifetime
            cap_repl_rate = -1 / lifetime

            # subtract the capital replacement (which is negative) rate
            # to the changes market share
            market_shares.loc[dict(region=region, year=self.year)] += cap_repl_rate

            # market decreasing faster than the average capital renewal rate
            # in this case, the idea is that oldest/non-competitive technologies
            # are likely to supply by increasing their lifetime
            # as the market does not justify additional capacity installation
            if volume_change < avg_cap_repl_rate:

                # we remove suppliers with a positive growth
                market_shares.loc[dict(region=region)].values[
                    market_shares.loc[dict(region=region)].values > 0
                ] = 0
                # we reverse the sign of negative growth suppliers
                market_shares.loc[dict(region=region)] *= -1
                market_shares.loc[dict(region=region)] /= market_shares.loc[
                    dict(region=region)
                ].sum(dim="variables")

                # multiply by volumes at T0
                market_shares.loc[dict(region=region)] *= data.sel(
                    region=region, year=self.year
                )
                market_shares.loc[dict(region=region)] /= market_shares.loc[
                    dict(region=region)
                ].sum(dim="variables")

            # increasing market or
            # market decreasing slowlier than the
            # capital renewal rate
            else:

                # we remove suppliers with a negative growth
                market_shares.loc[dict(region=region)].values[
                    market_shares.loc[dict(region=region)].values < 0
                ] = 0
                market_shares.loc[dict(region=region)] /= market_shares.loc[
                    dict(region=region)
                ].sum(dim="variables")

                # multiply by volumes at T0
                market_shares.loc[dict(region=region)] *= data.sel(
                    region=region, year=self.year
                )
                market_shares.loc[dict(region=region)] /= market_shares.loc[
                    dict(region=region)
                ].sum(dim="variables")

        return market_shares

    def __get_other_iam_vars(self, data: xr.DataArray) -> xr.DataArray:
        """
        Returns various IAM variables.

        :return: array containing various IAM variables
        """

        labels = self.__get_iam_variable_labels(IAM_OTHER_VARS, variable="iam_aliases")

        list_vars = list(labels.values())
        list_vars = [l for p in list_vars for l in p]

        data_to_return = data.loc[:, list_vars, :]

        return data_to_return

    def __get_iam_electricity_markets(self, data: xr.DataArray) -> xr.DataArray:
        """
        This method retrieves the market share for each electricity-producing technology,
        or a specified year, for each region provided by the IAM.

        :return: a multi-dimensional array with electricity technologies market share
        for a given year, for all regions.

        """

        labels = self.__get_iam_variable_labels(IAM_ELEC_VARS, variable="iam_aliases")

        list_technologies = list(labels.values())

        # If the year specified is not contained within the range of years given by the IAM
        if self.year < data.year.values.min() or self.year > data.year.values.max():
            raise KeyError(
                f"{self.year} is outside of the boundaries "
                f"of the IAM file: {data.year.values.min()}-{data.year.values.max()}"
            )

        # Finally, if the specified year falls in between two periods provided by the IAM
        # Interpolation between two periods

        try:
            data_to_return = data.loc[:, list_technologies, :]
        except KeyError as exc:
            list_missing_vars = [
                var for var in list_technologies if var not in data.variables.values
            ]
            print(
                f"The following variables cannot be found in the IAM file: {list_missing_vars}"
            )
            if len(list_technologies) - len(list_missing_vars) > 0:
                available_vars = [
                    var for var in list_technologies if var in data.variables.values
                ]
                print(
                    "The process continues with the remaining variables, "
                    "but certain transformation functions may not work."
                )
                list_technologies = available_vars
                data_to_return = data.loc[:, list_technologies, :]
            else:
                raise SystemExit from exc

        # give the array premise labels
        list_vars = [k for k, v in labels.items() if v in list_technologies]

        data_to_return.coords["variables"] = list_vars

        if self.system_model == "consequential":

            data_to_return = self.__transform_to_marginal_markets(data_to_return)

        else:
            data_to_return /= (
                data.loc[:, list_technologies, :].groupby("region").sum(dim="variables")
            )

        return data_to_return

    def __get_iam_electricity_efficiencies(self, data: xr.DataArray) -> xr.DataArray:
        """
        This method retrieves efficiency values for
        electricity-producing technology,
        for a specified year, for each region
        provided by the IAM.
        Electricity production from hydrogen can
        be removed from the mix
        (unless specified, it is removed).

        :return: a multi-dimensional array with electricity
        technologies market share for a given year, for all regions.

        """

        labels = self.__get_iam_variable_labels(IAM_ELEC_VARS, variable="eff_aliases")

        list_technologies = list(labels.values())

        # If the year specified is not contained within the range of years given by the IAM
        if self.year < data.year.values.min() or self.year > data.year.values.max():
            raise KeyError(
                f"{self.year} is outside of the boundaries "
                f"of the IAM file: {data.year.values.min()}-{data.year.values.max()}"
            )

        # Finally, if the specified year falls in between two periods provided by the IAM
        # Interpolation between two periods

        try:
            data_to_return = data.loc[:, list_technologies, :]
        except KeyError as exc:
            list_missing_vars = [
                var for var in list_technologies if var not in data.variables.values
            ]
            print(
                f"The following variables cannot be found in the IAM file: {list_missing_vars}"
            )
            if len(list_technologies) - len(list_missing_vars) > 0:
                available_vars = [
                    var for var in list_technologies if var in data.variables.values
                ]
                print(
                    "The process continues with the remaining variables, "
                    "but certain transformation functions may not work."
                )
                list_technologies = available_vars
                data_to_return = data.loc[:, list_technologies, :]
            else:
                raise SystemExit from exc

        data_to_return = data_to_return / data_to_return.sel(year=2020)

        # If we are looking at a year post 2020
        # and the ratio in efficiency change is inferior to 1
        # we correct it to 1, as we do not accept
        # that efficiency degrades over time
        data_to_return.loc[
            dict(year=[y for y in data_to_return.year.values if y > 2020])
        ] = np.clip(
            data_to_return.loc[
                dict(year=[y for y in data_to_return.year.values if y > 2020])
            ],
            1,
            None,
        )

        # Inversely, if we are looking at a year prior to 2020
        # and the ratio in efficiency change is superior to 1
        # we correct it to 1, as we do not accept
        # that efficiency in the past was higher than now
        data_to_return.loc[
            dict(year=[y for y in data_to_return.year.values if y < 2020])
        ] = np.clip(
            data_to_return.loc[
                dict(year=[y for y in data_to_return.year.values if y < 2020])
            ],
            None,
            1,
        )

        # ensure that efficiency does not decrease over time
        # ensure that efficiency can not decrease over time
        while data_to_return.diff(dim="year").min().values < 0:
            diff = data_to_return.diff(dim="year")
            diff = xr.concat([data_to_return.sel(year=2005), diff], dim="year")
            diff.values[diff.values > 0] = 0
            diff *= -1
            data_to_return += diff

        # convert NaNs to ones
        data_to_return = data_to_return.fillna(1)

        data_to_return.coords["variables"] = [
            k for k, v in labels.items() if v in list_technologies
        ]

        return data_to_return

    def __get_iam_cement_efficiencies(self, data: xr.DataArray) -> xr.DataArray:
        """
        This method retrieves specific energy use values for cement-producing technology,
        for a specified year, for each region provided by the IAM.

        :return: a multi-dimensional array with electricity technologies market share
        for a given year, for all regions.

        """

        # If the year specified is not contained within the range of years given by the IAM
        if self.year < data.year.values.min() or self.year > data.year.values.max():
            raise KeyError(
                f"{self.year} is outside of the boundaries "
                f"of the IAM file: {data.year.values.min()}-{data.year.values.max()}"
            )

        if (
            len(self.__get_iam_variable_labels(IAM_CEMENT_VARS, variable="eff_aliases"))
            > 0
        ):
            eff = self.__get_iam_variable_labels(
                IAM_CEMENT_VARS, variable="eff_aliases"
            )

            if eff["cement"] in data.variables.values:
                data_to_return = 1 / data.loc[:, [eff["cement"]], :]
            else:
                print("No efficiency variables is given for the cement sector.")
                data_to_return = xr.ones_like(data)
                var = data_to_return.variables.values.tolist()
                data_to_return = data_to_return.sel(variables=[var[0]])
        else:
            prod = self.__get_iam_variable_labels(
                IAM_CEMENT_VARS, variable="iam_aliases"
            )
            energy = self.__get_iam_variable_labels(
                IAM_CEMENT_VARS, variable="energy_use_aliases"
            )

            if (
                all(v in data.variables.values for v in energy["cement"])
                and prod["cement"] in data.variables.values
            ):

                data_to_return = 1 / (
                    data.loc[:, energy["cement"], :].sum(dim="variables")
                    / data.loc[:, [prod["cement"]], :]
                )
            else:
                print("No efficiency variables is given for the cement sector.")
                data_to_return = xr.ones_like(data)
                var = data_to_return.variables.values.tolist()
                data_to_return = data_to_return.sel(variables=[var[0]])

        data_to_return = data_to_return / data_to_return.sel(year=2020)

        # If we are looking at a year post 2020
        # and the ratio in specific energy use change is superior to 1
        # we correct it to 1, as we do not accept
        # that efficiency degrades over time
        data_to_return.loc[
            dict(year=[y for y in data_to_return.year.values if y > 2020])
        ] = np.clip(
            data_to_return.loc[
                dict(year=[y for y in data_to_return.year.values if y > 2020])
            ],
            1,
            None,
        )

        # Inversely, if we are looking at a year prior to 2020
        # and the ratio in specific energy use change is inferior to 1
        # we correct it to 1, as we do not accept
        # that efficiency in the past was higher than now
        data_to_return.loc[
            dict(year=[y for y in data_to_return.year.values if y < 2020])
        ] = np.clip(
            data_to_return.loc[
                dict(year=[y for y in data_to_return.year.values if y < 2020])
            ],
            None,
            1,
        )

        # ensure that efficiency does not decrease over time
        # ensure that efficiency can not decrease over time
        while data_to_return.diff(dim="year").min().values < 0:
            diff = data_to_return.diff(dim="year")
            diff = xr.concat([data_to_return.sel(year=2005), diff], dim="year")
            diff.values[diff.values > 0] = 0
            diff *= -1
            data_to_return += diff

        # convert NaNs to ones
        data_to_return = data_to_return.fillna(1)

        # we also consider any improvement rate
        # above 2 (+100%) or below 0.5 (-100%)
        # to be incorrect
        data_to_return.values = np.clip(data_to_return, 0.5, None)

        data_to_return.coords["variables"] = ["cement"]

        return data_to_return

    def __get_iam_steel_efficiencies(self, data: xr.DataArray) -> xr.DataArray:
        """
        This method retrieves specific energy use values for steel-producing technology,
        for a specified year, for each region provided by the IAM.

        :return: a multi-dimensional array with electricity technologies market share
        for a given year, for all regions.

        """

        # If the year specified is not contained within the range of years given by the IAM
        if self.year < data.year.values.min() or self.year > data.year.values.max():
            raise KeyError(
                f"{self.year} is outside of the boundaries "
                f"of the IAM file: {data.year.values.min()}-{data.year.values.max()}"
            )

        if (
            len(self.__get_iam_variable_labels(IAM_STEEL_VARS, variable="eff_aliases"))
            > 0
        ):
            eff = self.__get_iam_variable_labels(IAM_STEEL_VARS, variable="eff_aliases")

            if eff["steel - primary"] in data.variables.values:
                data_primary = 1 / data.loc[:, [eff["steel - primary"]], :]
            else:
                print("No efficiency variables is given for the primary steel sector.")
                data_primary = xr.ones_like(data)
                var = data_primary.variables.values.tolist()
                data_primary = data_primary.sel(variables=[var[0]])

        else:
            prod = self.__get_iam_variable_labels(
                IAM_STEEL_VARS, variable="iam_aliases"
            )
            energy = self.__get_iam_variable_labels(
                IAM_STEEL_VARS, variable="energy_use_aliases"
            )

            if isinstance(energy["steel - primary"], str):
                energy_in = [energy["steel - primary"]]
            else:
                energy_in = energy["steel - primary"]

            if prod["steel - primary"] in data.variables.values:
                data_primary = 1 / (
                    data.loc[:, energy_in, :].sum(dim="variables")
                    / data.loc[:, [prod["steel - primary"]], :]
                )
            else:
                print("No efficiency variables is given for the primary steel sector.")
                data_primary = xr.ones_like(data)
                var = data_primary.variables.values.tolist()
                data_primary = data_primary.sel(variables=[var[0]])

        # primary steel efficiency changes relative to 2020
        data_primary = data_primary / data_primary.sel(year=2020)

        if (
            len(self.__get_iam_variable_labels(IAM_STEEL_VARS, variable="eff_aliases"))
            > 0
        ):
            eff = self.__get_iam_variable_labels(IAM_STEEL_VARS, variable="eff_aliases")

            if eff["steel - secondary"] in data.variables.values:
                data_secondary = 1 / data.loc[:, [eff["steel - secondary"]], :]
            else:
                print(
                    "No efficiency variables is given for the secondary steel sector."
                )
                data_secondary = xr.ones_like(data)
                var = data_secondary.variables.values.tolist()
                data_secondary = data_secondary.sel(variables=[var[0]])

        else:
            prod = self.__get_iam_variable_labels(
                IAM_STEEL_VARS, variable="iam_aliases"
            )
            energy = self.__get_iam_variable_labels(
                IAM_STEEL_VARS, variable="energy_use_aliases"
            )

            if isinstance(energy["steel - secondary"], str):
                energy_in = [energy["steel - secondary"]]
            else:
                energy_in = energy["steel - secondary"]

            if prod["steel - secondary"] in data.variables.values:
                data_secondary = 1 / (
                    data.loc[:, energy_in, :].sum(dim="variables")
                    / data.loc[:, [prod["steel - secondary"]], :]
                )
            else:
                print(
                    "No efficiency variables is given for the secondary steel sector."
                )
                data_secondary = xr.ones_like(data)
                var = data_secondary.variables.values.tolist()
                data_secondary = data_secondary.sel(variables=[var[0]])

        # secondary steel efficiency changes relative to 2020
        data_secondary = data_secondary / data_secondary.sel(year=2020)

        data_to_return = xr.concat([data_primary, data_secondary], dim="variables")

        # If we are looking at a year post 2020
        # and the ratio in specific energy use change is superior to 1
        # we correct it to 1, as we do not accept
        # that efficiency degrades over time
        data_to_return.loc[
            dict(year=[y for y in data_to_return.year.values if y > 2020])
        ] = np.clip(
            data_to_return.loc[
                dict(year=[y for y in data_to_return.year.values if y > 2020])
            ],
            1,
            None,
        )

        # Inversely, if we are looking at a year prior to 2020
        # and the ratio in specific energy use change is inferior to 1
        # we correct it to 1, as we do not accept
        # that efficiency in the past was higher than now
        data_to_return.loc[
            dict(year=[y for y in data_to_return.year.values if y < 2020])
        ] = np.clip(
            data_to_return.loc[
                dict(year=[y for y in data_to_return.year.values if y < 2020])
            ],
            None,
            1,
        )

        # ensure that efficiency can not decrease over time
        while data_to_return.diff(dim="year").min().values < 0:
            diff = data_to_return.diff(dim="year")
            diff = xr.concat([data_to_return.sel(year=2005), diff], dim="year")
            diff.values[diff.values > 0] = 0
            diff *= -1
            data_to_return += diff

        # convert NaNs to ones
        data_to_return = data_to_return.fillna(1)

        # we also consider any improvement rate
        # above 2 (+100%) or below 0.5 (-100%)
        # to be incorrect
        data_to_return.values = np.clip(data_to_return, 0.5, None)

        data_to_return.coords["variables"] = [
            "steel - primary",
            "steel - secondary",
        ]

        return data_to_return

    def __get_gains_electricity_emissions(self, data: xr.DataArray) -> xr.DataArray:
        """
        This method retrieves emission values for electricity-producing technology,
        for a specified year, for each region provided by GAINS.

        :return: an multi-dimensional array with emissions for different technologies
        for a given year, for all regions.

        """

        labels = self.__get_iam_variable_labels(IAM_ELEC_VARS, variable="gains_aliases")

        # If the year specified is not contained within the range of years given by the IAM
        if self.year < data.year.values.min() or self.year > data.year.values.max():
            raise KeyError(
                f"{self.year} is outside of the boundaries "
                f"of the IAM file: {data.year.values.min()}-{data.year.values.max()}"
            )

        # Finally, if the specified year falls in between two periods provided by the IAM
        # Interpolation between two periods
        data_to_return = data.sel(sector=list(labels.values()))

        # Example: 5g CO per kWh in 2030, against 10g in 2020
        # 5/10 = 0.5
        # 1/0.5 = 2. Improvement factor of 2.
        data_to_return = 1 / (
            data_to_return.interp(year=self.year) / data_to_return.sel(year=2020)
        )

        # If we are looking at a year post 2020
        # and the ratio in efficiency change is inferior to 1
        # we correct it to 1, as we do not accept
        # that efficiency degrades over time
        if self.year > 2020:
            data_to_return.values[data_to_return.values < 1] = 1

        # Inversely, if we are looking at a year prior to 2020
        # and the ratio in efficiency change is superior to 1
        # we correct it to 1, as we do not accept
        # that efficiency in the past was higher than now
        if self.year < 2020:
            data_to_return.values[data_to_return.values > 1] = 1

        # convert NaNs to ones
        data_to_return = data_to_return.fillna(1)

        data_to_return.coords["sector"] = list(labels.keys())

        return data_to_return

    def __get_gains_cement_emissions(self, data: xr.DataArray) -> xr.DataArray:
        """
        This method retrieves emission values for cement production,
        for a specified year, for each region provided by GAINS.

        :return: a multi-dimensional array with emissions for different technologies
        for a given year, for all regions.


        """
        # If the year specified is not contained within the range of years given by the IAM
        if self.year < data.year.values.min() or self.year > data.year.values.max():
            raise KeyError(
                f"{self.year} is outside of the boundaries "
                f"of the IAM file: {data.year.values.min()}-{data.year.values.max()}"
            )

        # Finally, if the specified year falls in between two periods provided by the IAM
        # Interpolation between two periods
        data_to_return = data.sel(sector=["CEMENT"])

        # Example: 5g CO per kg cement in 2030, against 10g in 2020
        # 5/10 = 0.5
        # 1/0.5 = 2. Improvement factor of 2.
        data_to_return = 1 / (
            data_to_return.interp(year=self.year) / data_to_return.sel(year=2020)
        )

        # If we are looking at a year post 2020
        # and the ratio in efficiency change is inferior to 1
        # we correct it to 1, as we do not accept
        # that efficiency degrades over time
        if self.year > 2020:
            data_to_return.values[data_to_return.values < 1] = 1

        # Inversely, if we are looking at a year prior to 2020
        # and the ratio in efficiency change is superior to 1
        # we correct it to 1, as we do not accept
        # that efficiency in the past was higher than now
        if self.year < 2020:
            data_to_return.values[data_to_return.values > 1] = 1

        # convert NaNs to ones
        data_to_return = data_to_return.fillna(1)

        data_to_return.coords["sector"] = ["cement"]

        return data_to_return

    def __get_gains_steel_emissions(self, data: xr.DataArray) -> xr.DataArray:
        """
        This method retrieves emission values for steel production, for a specified year,
        for each region provided by GAINS.

        :return: a multi-dimensional array with emissions for different technologies
        for a given year, for all regions.

        """
        # If the year specified is not contained within the range of years given by the IAM
        if self.year < data.year.values.min() or self.year > data.year.values.max():
            raise KeyError(
                f"{self.year} is outside of the boundaries "
                f"of the IAM file: {data.year.values.min()}-{data.year.values.max()}"
            )

        # Finally, if the specified year falls in between two periods provided by the IAM
        # Interpolation between two periods
        data_to_return = data.sel(sector=["STEEL"])

        # Example: 5g CO per kg cement in 2030, against 10g in 2020
        # 5/10 = 0.5
        # 1/0.5 = 2. Improvement factor of 2.
        data_to_return = 1 / (
            data_to_return.interp(year=self.year) / data_to_return.sel(year=2020)
        )

        # If we are looking at a year post 2020
        # and the ratio in efficiency change is inferior to 1
        # we correct it to 1, as we do not accept
        # that efficiency degrades over time
        if self.year > 2020:
            data_to_return.values[data_to_return.values < 1] = 1

        # Inversely, if we are looking at a year prior to 2020
        # and the ratio in efficiency change is superior to 1
        # we correct it to 1, as we do not accept
        # that efficiency in the past was higher than now
        if self.year < 2020:
            data_to_return.values[data_to_return.values > 1] = 1

        # convert NaNs to ones
        data_to_return = data_to_return.fillna(1)

        data_to_return.coords["sector"] = ["steel"]

        return data_to_return

    def __get_iam_fuel_markets(self, data: xr.DataArray) -> xr.DataArray:
        """
        This method retrieves the market share
        for each fuel-producing technology,
        for a specified year, for each
        region provided by the IAM.

        :return: a multi-dimensional array with
        electricity technologies market share
        for a given year, for all regions.

        """

        labels = self.__get_iam_variable_labels(IAM_FUELS_VARS, variable="iam_aliases")

        list_technologies = list(labels.values())

        # If the year specified is not contained within the range of years given by the IAM
        if self.year < data.year.values.min() or self.year > data.year.values.max():
            raise KeyError(
                f"{self.year} is outside of the boundaries "
                f"of the IAM file: {data.year.values.min()}-{data.year.values.max()}"
            )

        # Finally, if the specified year falls in between two periods provided by the IAM
        # sometimes, the World region is either neglected
        # or wrongly evaluated, so we fix that here

        try:
            data.loc[dict(region="World", variables=list_technologies)] = data.loc[
                dict(
                    region=[r for r in data.coords["region"].values if r != "World"],
                    variables=list_technologies,
                )
            ].sum(dim="region")

        except KeyError as exc:
            list_missing_vars = [
                var for var in list_technologies if var not in data.variables.values
            ]
            print(
                f"The following variables cannot be found in the IAM file: {list_missing_vars}"
            )
            if len(list_technologies) - len(list_missing_vars) > 0:
                available_vars = [
                    var for var in list_technologies if var in data.variables.values
                ]
                print(
                    "The process continues with the remaining variables, "
                    "but certain transformation functions may not work."
                )
                list_technologies = available_vars
            else:
                raise SystemExit from exc

        # Interpolation between two periods
        data_to_return = data.loc[:, list_technologies, :]

        data_to_return.coords["variables"] = [
            k for k, v in labels.items() if v in list_technologies
        ]

        if self.system_model == "consequential":

            data_to_return = self.__transform_to_marginal_markets(data_to_return)

        else:
            data_to_return = data_to_return.interp(year=self.year)
            data_to_return /= (
                data.loc[:, list_technologies, :]
                .interp(year=self.year)
                .groupby("region")
                .sum(dim="variables")
            )

        return data_to_return

    def __get_iam_land_use(self, data):
        """
        Only provided by IMAGE at the moment. Those are land footprint
        associated with growing a given crop type, in hectares per GJ of that crop,
        for each region and year. This land occupation is added to the LCI
        for crop farming in fuels.py.

        :param data: IAM data
        :return: a multi-dimensional array with land use
        for different crops types, for all years, for all regions.
        """

        crops_vars = get_crops_properties()
        labels = list(crops_vars.keys())
        list_vars = [x["land_use"][self.model] for x in crops_vars.values()]

        data_to_return = data.loc[:, list_vars, :]
        data_to_return.coords["variables"] = labels

        return data_to_return

    def __get_iam_land_use_change_emissions(self, data: xr.DataArray) -> xr.DataArray:
        """
        Only provided by IMAGE at the moment. Those are CO2-eq. emissions
        associated with growing a given crop type, per GJ of that crop,
        for each region and year. Such LUC emissions are added to the LCI
        for crop farming in fuels.py.

        :param data: IAM data
        :return: a multi-dimensional array with land use change CO2 emissions
        for different crops types, for all years, for all regions.
        """

        crops_vars = get_crops_properties()
        labels = list(crops_vars.keys())
        list_vars = [x["land_use_change"][self.model] for x in crops_vars.values()]

        try:
            data_to_return = data.loc[:, list_vars, :]

        except KeyError:
            list_missing_vars = [
                var for var in list_vars if var not in data.variables.values
            ]
            raise KeyError(
                f"The following variables cannot be found in the IAM file: {list_missing_vars}"
            )

        data_to_return.coords["variables"] = labels

        return data_to_return

    def __get_iam_fuel_efficiencies(self, data: xr.DataArray) -> xr.DataArray:
        """
        This method retrieves the change in fuel production efficiency
        between the year in question and 2020, for each region provided by the IAM.
        Because we assume that the fuel conversion efficiency in ecoinvent or imported
        inventories are current (hence, representative of 2020).
        If the efficiency drops after 2020, we ignore it and keep the change
        in efficiency ratio to 1.

        :return: a multi-dimensional array with electricity technologies market
        share for a given year, for all regions.
        """

        labels = self.__get_iam_variable_labels(IAM_FUELS_VARS, variable="eff_aliases")

        list_technologies = list(labels.values())

        # If the year specified is not contained within the range of years given by the IAM
        if self.year < data.year.values.min() or self.year > data.year.values.max():
            raise KeyError(
                f"{self.year} is outside of the boundaries "
                f"of the IAM file: {data.year.values.min()}-{data.year.values.max()}"
            )

        # Finally, if the specified year falls in between two periods provided by the IAM
        # Interpolation between two periods
        try:
            data_to_return = data.loc[:, list_technologies, :]
        except KeyError as exc:
            list_missing_vars = [
                var for var in list_technologies if var not in data.variables.values
            ]
            print(
                f"The following variables cannot be found in the IAM file: {list_missing_vars}"
            )
            if len(list_technologies) - len(list_missing_vars) > 0:
                available_vars = [
                    var for var in list_technologies if var in data.variables.values
                ]
                print(
                    "The process continues with the remaining variables, "
                    "but certain transformation functions may not work."
                )
                list_technologies = available_vars
                data_to_return = data.loc[:, list_technologies, :]
            else:
                raise SystemExit from exc

        data_to_return = data_to_return / data_to_return.sel(year=2020)

        # If we are looking at a year post 2020
        # and the ratio in specific energy use change is superior to 1
        # we correct it to 1, as we do not accept
        # that efficiency degrades over time
        data_to_return.loc[
            dict(year=[y for y in data_to_return.year.values if y > 2020])
        ] = np.clip(
            data_to_return.loc[
                dict(year=[y for y in data_to_return.year.values if y > 2020])
            ],
            1,
            None,
        )

        # Inversely, if we are looking at a year prior to 2020
        # and the ratio in specific energy use change is inferior to 1
        # we correct it to 1, as we do not accept
        # that efficiency in the past was higher than now
        data_to_return.loc[
            dict(year=[y for y in data_to_return.year.values if y < 2020])
        ] = np.clip(
            data_to_return.loc[
                dict(year=[y for y in data_to_return.year.values if y < 2020])
            ],
            None,
            1,
        )

        # ensure that efficiency does not decrease over time
        # ensure that efficiency can not decrease over time
        while data_to_return.diff(dim="year").min().values < 0:
            diff = data_to_return.diff(dim="year")
            diff = xr.concat([data_to_return.sel(year=2005), diff], dim="year")
            diff.values[diff.values > 0] = 0
            diff *= -1
            data_to_return += diff

        # convert NaNs to ones
        data_to_return = data_to_return.fillna(1)

        # we also consider any improvement rate
        # above 2 (+100%) or below 0.5 (-100%)
        # to be incorrect
        data_to_return.values = np.clip(data_to_return.values, 0.5, None)

        data_to_return.coords["variables"] = [
            k for k, v in labels.items() if v in list_technologies
        ]

        return data_to_return

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

        cement_rate = data.loc[:, dict_vars["cement - cco2"], :].sum(
            dim=["variables"]
        ) / data.loc[:, dict_vars["cement - co2"], :].sum(dim=["variables"])

        cement_rate.coords["variables"] = "cement"

        steel_rate = data.loc[:, dict_vars["steel - cco2"], :].sum(
            dim="variables"
        ) / data.loc[:, dict_vars["steel - co2"], :].sum(dim="variables")
        steel_rate.coords["variables"] = "steel"

        rate = xr.concat([cement_rate, steel_rate], dim="variables")

        rate = rate.fillna(0)

        # we need to fix the rate for "World"
        # as it is sometimes neglected in the
        # IAM files

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

        # we ensure that the rate can only be between 0 and 1
        rate.values = np.clip(rate, 0, 1)

        return rate

    def __get_iam_production_volumes(self, dict_products, data) -> xr.DataArray:
        """
        Returns n xarray with production volumes for different sectors:
        electricity, steel, cement, fuels.
        This is used to build markets: we use
        the production volumes of each region for example,
        to build the World market.
        :param dict_products: a dictionary that contains
        common labels as keys, and IAM labels as values.
        :param data: IAM data
        :return: a xarray with production volumes for
        different commodities (electricity, cement, etc.)
        """

        list_products = list(dict_products.values())

        # If the year specified is not contained within the range of years given by the IAM
        if self.year < data.year.values.min() or self.year > data.year.values.max():
            raise KeyError(
                f"{self.year} is outside of the boundaries "
                f"of the IAM file: {data.year.values.min()}-{data.year.values.max()}"
            )

        # Finally, if the specified year falls in between two periods provided by the IAM
        # Interpolation between two periods

        try:

            data_to_return = data.loc[:, list_products, :]
        except KeyError as exc:
            list_missing_vars = [
                var for var in list_products if var not in data.variables.values
            ]
            print(
                f"The following variables cannot be found in the IAM file: {list_missing_vars}"
            )
            if len(list_products) - len(list_missing_vars) > 0:
                available_vars = [
                    var for var in list_products if var in data.variables.values
                ]
                print(
                    "The process continues with the remaining variables, "
                    "but certain transformation functions may not work."
                )
                list_products = available_vars
                data_to_return = data.loc[:, list_products, :]
            else:
                raise SystemExit from exc

        data_to_return.coords["variables"] = [
            k for k, v in dict_products.items() if v in list_products
        ]

        return data_to_return

    def get_external_data(self, datapackages):

        data = {}

        for i, dp in enumerate(datapackages):

            data[i] = {}

            resource = dp.get_resource("scenario_data")
            scenario_data = resource.read()
            scenario_headers = resource.headers
            df = pd.DataFrame(scenario_data, columns=scenario_headers)

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
                    (df["model"] == self.model)
                    & (df["pathway"] == self.pathway)
                    & (df["scenario"] == self.external_scenarios[i])
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

                array.coords["year"] = [int(y) for y in array.coords["year"]]

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
                        (df["model"] == self.model)
                        & (df["pathway"] == self.pathway)
                        & (df["scenario"] == self.external_scenarios[i])
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
                    array.coords["year"] = [int(y) for y in array.coords["year"]]

                    ref_years = {}

                    if "production pathways" in config_file:
                        for v in config_file["production pathways"].values():
                            for e, f in v.items():
                                if e == "efficiency":
                                    for x in f:
                                        ref_years[x["variable"]] = x.get(
                                            "reference year", 2020
                                        )

                    if "markets" in config_file:
                        for market in config_file["markets"]:
                            for e, f in market.items():
                                if f == "efficiency":
                                    for x in f["efficiency"]:
                                        ref_years[x["variable"]] = x.get(
                                            "reference year", 2020
                                        )

                    if ref_years:
                        for v, y in ref_years.items():

                            array.loc[dict(variables=v)] = array.loc[
                                dict(variables=v)
                            ] / array.loc[dict(variables=v)].sel(year=int(y))

                    # convert NaNs to ones
                    array = array.fillna(1)

                    data[i]["efficiency"] = array

        return data
