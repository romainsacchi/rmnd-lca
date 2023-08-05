"""
Integrates projections regarding use of metals in the economy,
from the DLR study:
Simon Schlichenmaier, Tobias Naegler,
May material bottlenecks hamper the global energy transition towards the 1.5°C target?,
Energy Reports, Volume 8, 2022, Pages 14875-14887, ISSN 2352-4847,
https://doi.org/10.1016/j.egyr.2022.11.025.
"""

import logging.config
import uuid
from functools import lru_cache
from pathlib import Path

import country_converter as coco
import numpy as np
import pandas as pd
import wurst
import yaml

from .export import biosphere_flows_dictionary
from .logger import create_logger
from .transformation import (
    BaseTransformation,
    Dict,
    IAMDataCollection,
    InventorySet,
    List,
    Set,
    ws,
)
from .utils import DATA_DIR

logger = create_logger("metal")


def _update_metals(scenario, version, system_model, modified_datasets):
    metals = Metals(
        database=scenario["database"],
        model=scenario["model"],
        pathway=scenario["pathway"],
        iam_data=scenario["iam data"],
        year=scenario["year"],
        version=version,
        system_model=system_model,
        modified_datasets=modified_datasets,
    )

    metals.create_metal_markets()

    return scenario, modified_datasets


def load_BGS_mapping():
    """
    Load mapping between BGS and ecoinvent
    """

    filepath = DATA_DIR / "metals" / "BGS_mapping.xlsx"
    df = pd.read_excel(filepath, sheet_name="Sheet1")
    return df


def get_ecoinvent_metal_factors():
    """
    Load dataframe with ecoinvent factors for metals
    and convert to xarray
    """

    filepath = DATA_DIR / "metals" / "ecoinvent_factors.csv"
    df = pd.read_csv(filepath)

    # create column "activity" as a tuple
    # of the columns "name", "product" and "location".
    df["activity"] = list(zip(df["name"], df["product"], df["location"]))
    df = df.drop(columns=["name", "product", "location"])
    df = df.melt(id_vars=["activity"], var_name="metal", value_name="value")

    # create an xarray with dimensions activity, metal and year
    ds = df.groupby(["activity", "metal"]).sum()["value"].to_xarray()

    return ds


def load_post_allocation_correction_factors():
    """
    Load yaml file with post-allocation correction factors
    """

    filepath = DATA_DIR / "metals" / "post-allocation correction" / "corrections.yaml"
    with open(filepath, "r", encoding="utf-8") as stream:
        factors = yaml.safe_load(stream)
    return factors


def fetch_mapping(filepath: str) -> dict:
    """Returns a dictionary from a YML file"""

    with open(filepath, "r", encoding="utf-8") as stream:
        mapping = yaml.safe_load(stream)
    return mapping


def rev_metals_map(mapping: dict) -> dict:
    """Returns a reversed dictionary"""

    rev_mapping = {}
    for key, val in mapping.items():
        for v in val:
            rev_mapping[v] = key
    return rev_mapping


def load_conversion_factors():
    """
    Load dataframe with conversion factors for metals
    """

    filepath = DATA_DIR / "metals" / "conversion_factors.xlsx"
    df = pd.read_excel(filepath, sheet_name="Conversion factors")
    return df


class Metals(BaseTransformation):
    """
    Class that modifies emissions of hot pollutants
    according to GAINS projections.
    """

    def __init__(
        self,
        database: List[dict],
        iam_data: IAMDataCollection,
        model: str,
        pathway: str,
        year: int,
        version: str,
        system_model: str,
        modified_datasets: dict,
    ):
        super().__init__(
            database,
            iam_data,
            model,
            pathway,
            year,
            version,
            system_model,
            modified_datasets,
        )

        self.version = version
        self.metals = iam_data.metals

        mapping = InventorySet(self.database, self.version)
        self.activities_metals_map: Dict[
            str, Set
        ] = mapping.generate_activities_using_metals_map()
        self.rev_activities_metals_map: Dict[str, str] = rev_metals_map(
            self.activities_metals_map
        )
        # self.metals_map: Dict[str, Set] = mapping.generate_metals_map()
        # self.rev_metals_map: Dict[str, str] = rev_metals_map(self.metals_map)
        self.conversion_factors = load_conversion_factors()
        self.current_metal_use = get_ecoinvent_metal_factors()

        self.biosphere_flow_codes = biosphere_flows_dictionary(version=self.version)

    def update_metals_use_in_database(self):
        """
        Update the database with metals use factors.
        """

        print("Integrating metals use factors.")
        for ds in self.database:
            if ds["name"] in self.rev_activities_metals_map:
                origin_var = self.rev_activities_metals_map[ds["name"]]
                self.update_metal_use(ds, origin_var)

            self.write_log(ds, "updated")

    def update_metal_use(
        self,
        dataset: dict,
        technology: str,
    ) -> dict:
        """
        Update metal use based on DLR data.
        :param dataset: dataset to adjust metal use for
        :param technology: DLR variable name to look up
        :return: Does not return anything. Modified in place.
        """

        # get the list of metal factors available for this technology

        if technology not in self.metals.origin_var.values:
            print(f"Technology {technology} not found in DLR data.")
            return dataset

        data = self.metals.sel(origin_var=technology, variable="median").interp(
            year=self.year
        )
        metals = [
            m
            for m in self.metals.metal.values
            if not np.isnan(data.sel(metal=m).values)
        ]

        # Update biosphere exchanges according to DLR use factors

        for exc in ws.biosphere(
            dataset, ws.either(*[ws.equals("name", x) for x in self.rev_metals_map])
        ):
            print("Updating metal use for", dataset["name"], exc["name"])
            metal = self.rev_metals_map[exc["name"]]
            use_factor = data.sel(metal=metal).values

            # check if there is a conversion factor
            if dataset["name"] in self.conversion_factors["Activity"].tolist():
                use_factor *= self.conversion_factors.loc[
                    self.conversion_factors["Activity"] == dataset["name"],
                    "Conversion_factor",
                ].values[0]

            else:
                print(f"Conversion factor not found for {dataset['name']}.")

            # update the exchange amount
            if metal in self.current_metal_use.metal.values:
                ecoinvent_factor = self.current_metal_use.sel(
                    metal=metal,
                    activity=(
                        dataset["name"],
                        dataset["reference product"],
                        dataset["location"],
                    ),
                ).values
            else:
                ecoinvent_factor = 0

            exc["amount"] += use_factor - ecoinvent_factor

            if "log parameters" not in dataset:
                dataset["log parameters"] = {}

            if metal not in dataset["log parameters"]:
                dataset["log parameters"][f"{metal} old amount"] = ecoinvent_factor
                dataset["log parameters"][f"{metal} new amount"] = exc["amount"]

                # remove metal from metals list
            metals.remove(metal)

        # Add new biosphere exchanges for metals
        # not present in the original dataset
        for metal in metals:
            use_factor = data.sel(metal=metal).values
            # check if there is a conversion factor
            if dataset["name"] in self.conversion_factors["Activity"].tolist():
                use_factor *= self.conversion_factors.loc[
                    self.conversion_factors["Activity"] == dataset["name"],
                    "Conversion_factor",
                ].values[0]
            else:
                print(f"Conversion factor not found for {dataset['name']}.")

            if self.version != "3.9":
                exc_id = (
                    f"{metal}, in ground",
                    "natural resource",
                    "in ground",
                    "kilogram",
                )
            else:
                exc_id = (
                    f"{metal}",
                    "natural resource",
                    "in ground",
                    "kilogram",
                )

            if metal in self.current_metal_use.metal.values:
                if (
                    dataset["name"],
                    dataset["reference product"],
                    dataset["location"],
                ) in self.current_metal_use.activity.values.tolist():
                    ecoinvent_factor = self.current_metal_use.sel(
                        metal=metal,
                        activity=(
                            dataset["name"],
                            dataset["reference product"],
                            dataset["location"],
                        ),
                    ).values
                else:
                    ecoinvent_factor = 0
            else:
                ecoinvent_factor = 0

            exc = {
                "name": f"{metal}, in ground",
                "amount": use_factor - ecoinvent_factor,
                "input": ("biosphere3", self.biosphere_flow_codes[exc_id]),
                "type": "biosphere",
                "unit": "kilogram",
                "comment": (f"{ecoinvent_factor};{use_factor};{technology};{metal}"),
            }

            dataset["exchanges"].append(exc)

            if "log parameters" not in dataset:
                dataset["log parameters"] = {}

            dataset["log parameters"][f"{metal} old amount"] = ecoinvent_factor
            dataset["log parameters"][f"{metal} new amount"] = exc["amount"]

        return dataset

    def post_allocation_correction(self):
        """
        Correct for post-allocation in the database.
        """

        factors_list = load_post_allocation_correction_factors()

        for dataset in factors_list:
            ds = ws.get_one(
                self.database,
                ws.equals("name", dataset["name"]),
                ws.equals("reference product", dataset["reference product"]),
                ws.equals("location", dataset["location"]),
                ws.equals("unit", dataset["unit"]),
            )
            ds["exchanges"].append(
                {
                    "name": dataset["additional flow"]["name"],
                    "amount": dataset["additional flow"]["amount"],
                    "unit": dataset["additional flow"]["unit"],
                    "type": "biosphere",
                    "categories": tuple(
                        dataset["additional flow"]["categories"].split("::")
                    ),
                    "input": (
                        "biosphere3",
                        self.biosphere_flow_codes[
                            dataset["additional flow"]["name"],
                            dataset["additional flow"]["categories"].split("::")[0],
                            dataset["additional flow"]["categories"].split("::")[1],
                            dataset["additional flow"]["unit"],
                        ],
                    ),
                }
            )

    def create_new_mining_activity(
        self, name, reference_product, new_locations, geo_mapping
    ) -> dict:
        """
        Create a new mining activity in a new location.
        """
        # Get the original datasets
        datasets = self.fetch_proxies(
            name=name,
            ref_prod=reference_product,
            regions=new_locations.values(),
            geo_mapping=geo_mapping,
            exact_match=True,
        )

        return datasets

    @lru_cache
    def convert_long_to_short_country_name(self, country_long):
        """
        Convert long country name to short country name.
        :param country_long: Long country name
        :return: Short country name
        """
        country_short = coco.convert(country_long, to="ISO2")

        if isinstance(country_short, list):
            if country_long == "France (French Guiana)":
                country_short = "GF"
            else:
                print(f"Multiple locations found for {country_long}. Using first one.")
                country_short = country_short[0]

        if country_short not in self.geo.geo.keys():
            print(f"New location {country_short} for {country_long} not found")
            return None

        return country_short

    def get_shares(self, df: pd.DataFrame, new_locations: dict, name, ref_prod) -> dict:
        """
        Get shares of each location in the dataframe.
        :param df: Dataframe with mining shares
        :param new_locations: List of new locations
        :return: Dictionary with shares of each location
        """
        shares = {}

        for long_location, short_location in new_locations.items():
            shares[(name, ref_prod, short_location)] = df.loc[
                df["Country"] == long_location, "Share_2017_2021"
            ].values[0]

        return shares

    def get_geo_mapping(self, df: pd.DataFrame, new_locations: dict) -> dict:
        """
        Fetch the value under "Region" for each location in new_locations.
        """

        mapping = {}

        for long_location, short_location in new_locations.items():
            mapping[short_location] = df.loc[
                df["Country"] == long_location, "Region"
            ].values[0]

        return mapping

    def create_region_specific_markets(self, df: pd.DataFrame) -> List[dict]:
        new_exchanges = []

        # iterate through unique pair of process - reference product in df:
        # for each pair, create a new mining activity

        for (name, ref_prod), group in df.groupby(["Process", "Reference product"]):
            new_locations = {
                c: self.convert_long_to_short_country_name(c)
                for c in group["Country"].unique()
            }
            # remove None values
            new_locations = {k: v for k, v in new_locations.items() if v}
            # fetch shares for each location in df
            shares = self.get_shares(group, new_locations, name, ref_prod)

            geography_mapping = self.get_geo_mapping(group, new_locations)

            # if not, we create it
            datasets = self.create_new_mining_activity(
                name, ref_prod, new_locations, geography_mapping
            )

            # add new datasets to database
            self.database.extend(datasets.values())

            self.modified_datasets[(self.model, self.scenario, self.year)][
                "created"
            ].extend(
                [
                    (
                        dataset["name"],
                        dataset["reference product"],
                        dataset["location"],
                        dataset["unit"],
                    )
                    for dataset in datasets.values()
                ]
            )

            new_exchanges.extend(
                [
                    {
                        "name": dataset["name"],
                        "product": dataset["reference product"],
                        "location": dataset["location"],
                        "unit": dataset["unit"],
                        "amount": shares[(name, ref_prod, dataset["location"])],
                        "type": "technosphere",
                    }
                    for dataset in datasets.values()
                ]
            )

        return new_exchanges

    def create_market(self, metal, df):
        # check if market already exists
        # if so, remove it

        for ds in self.database:
            if ds["name"] == f"market for {metal[0].lower() + metal[1:]}":
                # add it to self.modified_datasets
                self.modified_datasets[(self.model, self.scenario, self.year)][
                    "emptied"
                ].append(
                    (
                        ds["name"],
                        ds["reference product"],
                        ds["location"],
                        ds["unit"],
                    )
                )
                # self.database.remove(ds)

        dataset = {
            "name": f"market for {metal[0].lower() + metal[1:]}",
            "location": "World",
            "exchanges": [
                {
                    "name": f"market for {metal[0].lower() + metal[1:]}",
                    "product": f"{metal[0].lower() + metal[1:]}",
                    "location": "World",
                    "amount": 1,
                    "type": "production",
                    "unit": "kilogram",
                }
            ],
            "reference product": f"{metal[0].lower() + metal[1:]}",
            "unit": "kilogram",
            "production amount": 1,
            "comment": "Created by premise",
            "code": str(uuid.uuid4()),
        }

        dataset["exchanges"].extend(self.create_region_specific_markets(df))

        # filter out None
        dataset["exchanges"] = [exc for exc in dataset["exchanges"] if exc]

        return dataset

    def create_metal_markets(self):
        self.post_allocation_correction()

        print("Creating metal markets")

        dataframe = load_BGS_mapping()
        dataframe = dataframe.loc[dataframe["Work done"] == "Yes"]
        dataframe = dataframe.loc[~dataframe["Country"].isnull()]
        dataframe_shares = dataframe.loc[dataframe["Share_2017_2021"] > 0]

        for metal in dataframe_shares["Metal"].unique():
            print(f"... for {metal}.")
            df_metal = dataframe.loc[dataframe["Metal"] == metal]
            dataset = self.create_market(metal, df_metal)

            self.database.append(dataset)

            # add it to self.modified_datasets
            self.modified_datasets[(self.model, self.scenario, self.year)][
                "created"
            ].append(
                (
                    dataset["name"],
                    dataset["reference product"],
                    dataset["location"],
                    dataset["unit"],
                )
            )

        dataframe_parent = dataframe.loc[
            (dataframe["Share_2017_2021"].isnull()) & (~dataframe["Region"].isnull())
        ]

        print("Creating additional mining processes")
        for metal in dataframe_parent["Metal"].unique():
            df_metal = dataframe_parent.loc[dataframe["Metal"] == metal]

            for (name, ref_prod), group in df_metal.groupby(
                ["Process", "Reference product"]
            ):
                new_locations = {
                    c: self.convert_long_to_short_country_name(c)
                    for c in group["Country"].unique()
                }
                # remove None
                new_locations = {
                    k: v for k, v in new_locations.items() if v is not None
                }

                geography_mapping = self.get_geo_mapping(group, new_locations)

                # if not, we create it
                datasets = self.create_new_mining_activity(
                    name, ref_prod, new_locations, geography_mapping
                )

                self.database.extend(datasets.values())

                # add it to self.modified_datasets
                self.modified_datasets[(self.model, self.scenario, self.year)][
                    "created"
                ].extend(
                    [
                        (
                            dataset["name"],
                            dataset["reference product"],
                            dataset["location"],
                            dataset["unit"],
                        )
                        for dataset in datasets.values()
                    ]
                )

    def write_log(self, dataset, status="created"):
        """
        Write log file.
        """

        if "log parameters" in dataset:
            logger.info(
                f"{status}|{self.model}|{self.scenario}|{self.year}|"
                f"{dataset['name']}|{dataset['location']}|"
            )
