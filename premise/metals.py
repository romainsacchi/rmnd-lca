"""
Integrates projections regarding use of metals in the economy,
from the DLR study:
Simon Schlichenmaier, Tobias Naegler,
May material bottlenecks hamper the global energy transition towards the 1.5°C target?,
Energy Reports, Volume 8, 2022, Pages 14875-14887, ISSN 2352-4847,
https://doi.org/10.1016/j.egyr.2022.11.025.
"""

from pathlib import Path
import logging.config

import numpy as np
import pandas as pd
import yaml

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

logger = logging.getLogger("emissions")


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
        super().__init__(database, iam_data, model, pathway, year, version, system_model, modified_datasets)

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

        self.biosphere_flow_codes = biosphere_flows_dictionary(
            version=self.version
        )

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


    def write_log(self, dataset, status="created"):
        """
        Write log file.
        """

        if "log parameters" in dataset:
            logger.info(
                f"{status}|{self.model}|{self.scenario}|{self.year}|"
                f"{dataset['name']}|{dataset['location']}|"
            )
