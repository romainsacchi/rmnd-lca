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
from pathlib import Path

import numpy as np
import pandas as pd
import wurst
import yaml
import country_converter as coco

from .export import biosphere_flows_dictionary
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

LOG_CONFIG = DATA_DIR / "utils" / "logging" / "logconfig.yaml"
# directory for log files
DIR_LOG_REPORT = Path.cwd() / "export" / "logs"
# if DIR_LOG_REPORT folder does not exist
# we create it
if not Path(DIR_LOG_REPORT).exists():
    Path(DIR_LOG_REPORT).mkdir(parents=True, exist_ok=True)

with open(LOG_CONFIG, "r") as f:
    config = yaml.safe_load(f.read())
    logging.config.dictConfig(config)

logger = logging.getLogger("metal")


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
        self.metals_map: Dict[str, Set] = mapping.generate_metals_map()
        self.rev_metals_map: Dict[str, str] = rev_metals_map(self.metals_map)
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
                    "categories": tuple(dataset["additional flow"]["categories"].split("::")),
                    "input": (
                        "biosphere3",
                        self.biosphere_flow_codes[
                            dataset["additional flow"]["name"],
                            dataset["additional flow"]["categories"].split("::")[0],
                            dataset["additional flow"]["categories"].split("::")[1],
                            dataset["additional flow"]["unit"],
                        ]
                    ),
                }
            )

    def create_new_mining_activity(self, name, reference_product, location, new_location):

        try:
            ws.get_one(
                self.database,
                ws.equals("name", name),
                ws.equals("reference product", reference_product),
                ws.equals("location", new_location),
            )

            return None

        except ws.NoResults:

            try:
                original = ws.get_one(
                    self.database,
                    ws.equals("name", name),
                    ws.equals("reference product", reference_product),
                    ws.equals("location", location),
                )
            except ws.NoResults:
                print(f"No original dataset found for {name}, {reference_product} in {location}")
                return None
            except ws.MultipleResults:
                print(f"Multiple original datasets found for {name}, {reference_product} in {location}")
                return None

            dataset = wurst.copy_to_new_location(original, new_location)
            dataset = self.relink_technosphere_exchanges(dataset)
            dataset["code"] = str(uuid.uuid4())

            return dataset

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

    def create_region_specific_market(self, row):

        new_location = self.convert_long_to_short_country_name(row["Country"])

        if new_location is None:
            return None

        name = row["Process"]
        reference_product = row["Reference product"]
        location = row["Region"]

        # check if already exists
        try:
            ds = ws.get_one(
                self.database,
                ws.equals("name", name),
                ws.equals("reference product", reference_product),
                ws.equals("location", new_location),
            )
            return {
                "name": ds["name"],
                "product": ds["reference product"],
                "location": ds["location"],
                "unit": ds["unit"],
                "amount": float(row["Share_2017_2021"]),
                "uncertainty type": 0,
                "type": "technosphere",
            }
        except ws.NoResults:
            pass

        # if not, we create it
        dataset = self.create_new_mining_activity(
            name,
            reference_product,
            location,
            new_location
        )

        if dataset is None:
            return None

        # add new dataset to database
        self.database.append(dataset)

        self.modified_datasets[
            (
                self.model,
                self.scenario,
                self.year
            )
        ][
            "created"
        ].append(
            (
                dataset["name"],
                dataset["reference product"],
                dataset["location"],
                dataset["unit"],
            )
        )

        return {
            "name": dataset["name"],
            "product": dataset["reference product"],
            "location": dataset["location"],
            "unit": dataset["unit"],
            "amount": float(row["Share_2017_2021"]),
            "uncertainty type": 0,
            "type": "technosphere",
        }

    def create_market(self, metal, df):

        # check if market already exists
        # if so, remove it
        for ds in self.database:
            if ds['name'] == f"market for {metal.lower()}":
                self.database.remove(ds)
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


        dataset = {
            "name": f"market for {metal.lower()}",
            "location": "World",
            "exchanges": [
                {
                    "name": f"market for {metal.lower()}",
                    "product": f"{metal.lower()}",
                    "location": "World",
                    "amount": 1,
                    "type": "production",
                    "unit": "kilogram",
                }
            ],
            "reference product": f"{metal.lower()}",
            "unit": "kilogram",
            "production amount": 1,
            "comment": "Created by premise",
            "code": str(uuid.uuid4()),
        }

        dataset["exchanges"].extend(
            [
                self.create_region_specific_market(row)
                for _, row in df.iterrows()
            ]
        )

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
            (dataframe["Share_2017_2021"].isnull())
            & (~dataframe["Region"].isnull())
        ]

        print("Creating additional mining processes")
        for _, row in dataframe_parent.iterrows():
            dataset = self.create_new_mining_activity(
                row["Process"],
                row["Reference product"],
                row["Region"],
                self.convert_long_to_short_country_name(row["Country"])
            )

            if dataset is None:
                continue

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


    def write_log(self, dataset, status="created"):
        """
        Write log file.
        """

        if "log parameters" in dataset:
            logger.info(
                f"{status}|{self.model}|{self.scenario}|{self.year}|"
                f"{dataset['name']}|{dataset['location']}|"
            )
