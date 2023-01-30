"""
Integrates projections regarding use of metals in the economy,
from the DLR study:
Simon Schlichenmaier, Tobias Naegler,
May material bottlenecks hamper the global energy transition towards the 1.5°C target?,
Energy Reports, Volume 8, 2022, Pages 14875-14887, ISSN 2352-4847,
https://doi.org/10.1016/j.egyr.2022.11.025.
"""

from datetime import date
from pathlib import Path

import pandas as pd
import numpy as np
import yaml

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
from .export import biosphere_flows_dictionary

EI_METALS = DATA_DIR / "metals" / "ecoinvent_metals.yaml"
LOG_DIR = DATA_DIR / "logs"

biosphere_flow_codes = biosphere_flows_dictionary()


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
    ):
        super().__init__(database, iam_data, model, pathway, year)

        self.version = version
        self.metals = iam_data.metals
        self.ei_metals = fetch_mapping(EI_METALS)

        mapping = InventorySet(self.database)
        self.metals_map: Dict[str, Set] = mapping.generate_metals_map()
        self.rev_metals_map: Dict[str, str] = rev_metals_map(self.metals_map)
        self.conversion_factors = load_conversion_factors()

    def update_metals_use_in_database(self):
        """
        Update the database with metals use factors.
        """

        print("Integrating metals use factors.")
        for ds in self.database:
            if ds["name"] in self.rev_metals_map:
                origin_var = self.rev_metals_map[ds["name"]]
                self.update_metal_use(ds, origin_var)

        self.logging_changes()

    def update_metal_use(
        self, dataset: dict, technology: str,
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

        data = self.metals.sel(origin_var=technology, variable="mean").interp(year=self.year)
        metals = [m for m in self.metals.metal.values if not np.isnan(data.sel(metal=m).values)]

        # Update biosphere exchanges according to DLR use factors
        for exc in ws.biosphere(
            dataset, ws.either(*[ws.equals("name", x) for x in self.ei_metals])
        ):

            metal = self.ei_metals[exc["name"]]
            use_factor = data.sel(metal=metal).values

            # check if there is a conversion factor
            if dataset["name"] in self.conversion_factors["Activity"].tolist():
                use_factor *= self.conversion_factors.loc[
                    self.conversion_factors["Activity"] == dataset["name"],
                    "Conversion_factor"
                ].values[0]
            else:
                print(
                    f"Conversion factor not found for {dataset['name']}."
                )

            if metal not in exc.get("comment", ""):
                exc["comment"] = (
                    f"{exc['amount']};{use_factor};{technology};{metal}"
                )

            # update the exchange amount
            exc["amount"] = use_factor
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
                    "Conversion_factor"
                ].values[0]
            else:
                print(
                    f"Conversion factor not found for {dataset['name']}."
                )
            exc_id = (
                f"{self.ei_metals[metal]}, in ground",
                "natural resource",
                "in ground",
                "kilogram"
            )

            exc = {
                "name": f"{self.ei_metals[metal]}, in ground",
                "amount": use_factor,
                "input": ("biosphere3", biosphere_flow_codes[exc_id]),
                "type": "biosphere",
                "unit": "kilogram",
                "comment": (
                    f"{use_factor};{use_factor};{technology};{metal}"
                )
            }
            dataset["exchanges"].append(exc)

        return dataset

    def logging_changes(self):
        """
        Log changes made to the database.
        """

        print("Logging changes made to the database.")
        list_res = []
        for ds in self.database:
            for exc in ds["exchanges"]:
                if exc["type"] == "biosphere" and exc.get("comment"):
                    if len(exc["comment"].split(";")) == 4:
                        d = [
                            ds["name"],
                            ds["reference product"],
                            ds["location"],
                            exc["name"],
                        ]
                        d.extend(exc["comment"].split(";"))
                        list_res.append(d)

        modified = pd.DataFrame(
            list_res,
            columns=[
                "dataset name",
                "dataset product",
                "dataset location",
                "metal",
                "old value",
                "new value",
                "DLR variable",
                "DLR metal"
            ],
        )

        # check that directory exists, otherwise create it
        Path(LOG_DIR).mkdir(parents=True, exist_ok=True)
        modified.to_csv(
            LOG_DIR
            / f"log modified metals use {self.model.upper()} {self.scenario} {self.year}-{date.today()}.csv",
            index=False,
        )
        print(f"Log of changes in metals use saved in {LOG_DIR}.")
