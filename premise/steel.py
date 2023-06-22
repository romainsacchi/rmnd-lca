"""
Integrates projections regarding steel production.
"""
import logging.config
from pathlib import Path
from typing import Dict, List

import wurst
import yaml

from .data_collection import IAMDataCollection
from .transformation import BaseTransformation, ws
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

logger = logging.getLogger("steel")


class Steel(BaseTransformation):
    """
    Class that modifies steel markets in ecoinvent based on IAM output data.

    :ivar database: database dictionary from :attr:`.NewDatabase.database`
    :ivar model: can be 'remind' or 'image'. str from :attr:`.NewDatabase.model`
    :ivar iam_data: xarray that contains IAM data, from :attr:`.NewDatabase.rdc`
    :ivar year: year, from :attr:`.NewDatabase.year`

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
    ) -> None:
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

    def generate_activities(self):
        """
        This function generates new activities for primary and secondary steel
        production and add them to the wurst database.

        :return: Returns a modified database with newly added steel activities for the corresponding year
        """

        self.create_steel_markets()
        self.create_steel_production_activities()
        self.create_pig_iron_markets()
        self.create_pig_iron_production_activities()

    def create_steel_markets(self):
        """
        Create steel markets for different regions

        :return: Does not return anything. Adds new markets to database.
        """

        steel_markets_to_create = (
            ("market for steel, low-alloyed", "steel, low-alloyed"),
            ("market for steel, unalloyed", "steel, unalloyed"),
            ("market for steel, chromium steel 18/8", "steel, chromium steel 18/8"),
        )

        list_new_steel_markets = []

        for market, steel_product in steel_markets_to_create:
            steel_markets = self.fetch_proxies(
                name=market,
                ref_prod=steel_product,
                production_variable=["steel - primary", "steel - secondary"],
            )

            # adjust share of primary and secondary steel
            if market == "market for steel, low-alloyed":
                for loc, dataset in steel_markets.items():
                    if loc != "World":
                        if self.system_model != "consequential":
                            try:
                                primary_share = self.iam_data.production_volumes.sel(
                                    region=loc, variables="steel - primary"
                                ).interp(year=self.year).values.item(
                                    0
                                ) / self.iam_data.production_volumes.sel(
                                    region=loc,
                                    variables=["steel - primary", "steel - secondary"],
                                ).interp(
                                    year=self.year
                                ).sum(
                                    dim="variables"
                                ).values.item(
                                    0
                                )
                            except KeyError:
                                primary_share = 1
                        else:
                            primary_share = 1

                        secondary_share = 1 - primary_share

                        new_exc = [
                            {
                                "uncertainty type": 0,
                                "loc": primary_share,
                                "amount": primary_share,
                                "type": "technosphere",
                                "production volume": 1,
                                "product": steel_product,
                                "name": "steel production, converter, low-alloyed",
                                "unit": "kilogram",
                                "location": loc,
                            }
                        ]

                        if secondary_share > 0:
                            new_exc.append(
                                {
                                    "uncertainty type": 0,
                                    "loc": secondary_share,
                                    "amount": secondary_share,
                                    "type": "technosphere",
                                    "production volume": 1,
                                    "product": steel_product,
                                    "name": "steel production, electric, low-alloyed",
                                    "unit": "kilogram",
                                    "location": loc,
                                },
                            )

                        dataset["exchanges"] = [
                            e
                            for e in dataset["exchanges"]
                            if e["type"] == "production" or e["unit"] == "ton kilometer"
                        ]
                        dataset["exchanges"].extend(new_exc)

                        dataset["log parameters"] = {
                            "primary steel share": primary_share,
                            "secondary steel share": secondary_share,
                        }
            else:
                for loc, dataset in steel_markets.items():
                    if loc != "World":
                        name_ref = [
                            (e["name"], e.get("product"))
                            for e in dataset["exchanges"]
                            if "steel production" in e["name"]
                        ][0]
                        name, ref = name_ref

                        dataset["exchanges"] = [
                            e
                            for e in dataset["exchanges"]
                            if e["type"] == "production" or e["unit"] == "ton kilometer"
                        ]

                        dataset["exchanges"].append(
                            {
                                "uncertainty type": 0,
                                "loc": 1,
                                "amount": 1,
                                "type": "technosphere",
                                "production volume": 1,
                                "product": ref,
                                "name": name,
                                "unit": "kilogram",
                                "location": loc,
                            }
                        )

            # populate World dataset
            steel_markets["World"]["exchanges"] = [
                x
                for x in steel_markets["World"]["exchanges"]
                if x["type"] == "production"
            ]
            regions = [r for r in self.regions if r != "World"]

            for region in regions:
                try:
                    share = (
                        self.iam_data.production_volumes.sel(
                            variables=["steel - primary", "steel - secondary"],
                            region=region,
                        )
                        .interp(year=self.year)
                        .sum(dim="variables")
                        / self.iam_data.production_volumes.sel(
                            variables=["steel - primary", "steel - secondary"],
                            region="World",
                        )
                        .interp(year=self.year)
                        .sum(dim="variables")
                    ).values.item(0)

                except KeyError:
                    # equal share to all regions
                    share = 1 / len(regions)

                steel_markets["World"]["exchanges"].append(
                    {
                        "name": market,
                        "product": steel_product,
                        "amount": share,
                        "unit": "kilogram",
                        "type": "technosphere",
                        "location": region,
                    }
                )

            # add to log
            for new_dataset in list(steel_markets.values()):
                self.write_log(new_dataset)
                # add it to list of created datasets
                self.modified_datasets[(self.model, self.scenario, self.year)][
                    "created"
                ].append(
                    (
                        new_dataset["name"],
                        new_dataset["reference product"],
                        new_dataset["location"],
                        new_dataset["unit"],
                    )
                )

            list_new_steel_markets.extend(list(steel_markets.values()))

        # add new steel markets to database
        self.database.extend(list_new_steel_markets)

    def create_steel_production_activities(self):
        """
        Create steel production activities for different regions.

        """
        # Determine all steel activities in the database. Empty old datasets.
        print("Create new steel production datasets and empty old datasets")

        d_act_primary_steel = {
            mat: self.fetch_proxies(
                name=mat[0],
                ref_prod=mat[1],
                production_variable=["steel - primary"],
            )
            for mat in zip(
                self.material_map["steel - primary"],
                ["steel"] * len(self.material_map["steel - primary"]),
            )
        }

        d_act_secondary_steel = {
            mat: self.fetch_proxies(
                name=mat[0],
                ref_prod=mat[1],
                production_variable=["steel - secondary"],
            )
            for mat in zip(
                self.material_map["steel - secondary"],
                ["steel"] * len(self.material_map["steel - secondary"]),
            )
        }

        # adjust efficiency of primary steel production
        # and add carbon capture and storage, if needed
        for _, steel in d_act_primary_steel.items():
            steel = self.adjust_process_efficiency(steel)
            steel = self.add_carbon_capture_and_storage(
                datasets=steel,
            )
            # update the database with the modified datasets
            self.database.extend(list(steel.values()))

            # add to log
            for new_dataset in list(steel.values()):
                self.write_log(new_dataset)
                # add it to list of created datasets
                self.modified_datasets[(self.model, self.scenario, self.year)][
                    "created"
                ].append(
                    (
                        new_dataset["name"],
                        new_dataset["reference product"],
                        new_dataset["location"],
                        new_dataset["unit"],
                    )
                )

        # adjust efficiency of secondary steel production
        # and add carbon capture and storage, if needed
        for _, steel in d_act_secondary_steel.items():
            steel = self.adjust_process_efficiency(steel)
            steel = self.add_carbon_capture_and_storage(
                datasets=steel,
            )
            # update the database with the modified datasets
            self.database.extend(list(steel.values()))

            # add to log
            for new_dataset in list(steel.values()):
                self.write_log(new_dataset)
                # add it to list of created datasets
                self.modified_datasets[(self.model, self.scenario, self.year)][
                    "created"
                ].append(
                    (
                        new_dataset["name"],
                        new_dataset["reference product"],
                        new_dataset["location"],
                        new_dataset["unit"],
                    )
                )

    def create_pig_iron_production_activities(self):
        """
        Create region-specific pig iron production activities.
        """

        print("Create pig iron production datasets")

        pig_iron = self.fetch_proxies(
            name="pig iron production",
            ref_prod="pig iron",
            production_variable=["steel - primary"],
        )

        # adjust efficiency of pig iron production
        pig_iron = self.adjust_process_efficiency(pig_iron)
        # add carbon capture and storage, if needed
        pig_iron = self.add_carbon_capture_and_storage(
            datasets=pig_iron,
        )

        self.database.extend(list(pig_iron.values()))

        # add to log
        for new_dataset in list(pig_iron.values()):
            self.write_log(new_dataset)
            # add it to list of created datasets
            self.modified_datasets[(self.model, self.scenario, self.year)][
                "created"
            ].append(
                (
                    new_dataset["name"],
                    new_dataset["reference product"],
                    new_dataset["location"],
                    new_dataset["unit"],
                )
            )

    def create_pig_iron_markets(self):
        """
        Create region-specific pig iron markets.
        Adds datasets to the database.
        """

        pig_iron_markets = self.fetch_proxies(
            name="market for pig iron",
            ref_prod="pig iron",
            production_variable=["steel - primary"],
        )
        self.database.extend(list(pig_iron_markets.values()))

        # add to log
        for new_dataset in list(pig_iron_markets.values()):
            self.write_log(new_dataset)
            # add it to list of created datasets
            self.modified_datasets[(self.model, self.scenario, self.year)][
                "created"
            ].append(
                (
                    new_dataset["name"],
                    new_dataset["reference product"],
                    new_dataset["location"],
                    new_dataset["unit"],
                )
            )

    def adjust_process_efficiency(self, datasets):
        """
        Scale down fuel exchanges in the given datasets, according to efficiency improvement as
        forecast by the IAM.

        :param datasets: A dictionary of datasets to modify.
        :type datasets: dict
        :return: The modified datasets.
        :rtype: dict
        """
        list_fuels = [
            "diesel",
            "coal",
            "lignite",
            "coke",
            "fuel",
            "meat",
            "gas",
            "oil",
            "electricity",
            "natural gas",
            "steam",
        ]

        for region, dataset in datasets.items():
            # Determine the sector based on the activity name
            if any(i in dataset["name"] for i in ["converter", "pig iron"]):
                sector = "steel - primary"
            else:
                sector = "steel - secondary"

            if sector in self.iam_data.steel_efficiencies.variables.values:
                # Calculate the scaling factor based on the efficiency change from 2020 to the current year
                scaling_factor = 1 / self.find_iam_efficiency_change(
                    data=self.iam_data.steel_efficiencies,
                    variable=sector,
                    location=dataset["location"],
                )
            else:
                scaling_factor = 1

            if scaling_factor != 1:
                # Update the comments
                text = (
                    f"This dataset has been modified by `premise`, according to the performance "
                    f"for steel production indicated by the IAM model {self.model.upper()} for the IAM "
                    f"region {region} in {self.year}, following the scenario {self.scenario}. "
                    f"The energy efficiency of the process has been improved by {int((1 - scaling_factor) * 100)}%."
                )
                dataset["comment"] = text + dataset["comment"]

                # Scale down the fuel exchanges using the scaling factor
                wurst.change_exchanges_by_constant_factor(
                    dataset,
                    scaling_factor,
                    technosphere_filters=[
                        ws.either(*[ws.contains("name", x) for x in list_fuels])
                    ],
                    biosphere_filters=[ws.contains("name", "Carbon dioxide, fossil")],
                )

                if "log parameters" not in dataset:
                    dataset["log parameters"] = {}

                dataset["log parameters"].update(
                    {
                        "thermal efficiency change": scaling_factor,
                    }
                )

        return datasets

    def add_carbon_capture_and_storage(self, datasets: Dict[str, dict]):
        """
        Adds carbon capture-related energy exchanges to the input datasets for the given sector.

        :param datasets: A dictionary of datasets to modify.
        :return: A modified dictionary of datasets.
        """

        for region, dataset in datasets.items():
            # Check if carbon capture rate data is available for this region and sector
            carbon_capture_rate = self.get_carbon_capture_rate(
                loc=dataset["location"], sector="steel"
            )
            if carbon_capture_rate > 0:
                # Create a new CCS dataset if one doesn't exist
                self.create_ccs_dataset(
                    loc=region, bio_co2_stored=0, bio_co2_leaked=0, sector="steel"
                )

                # Modify the CO2 flow in the input dataset
                for co2_flow in ws.biosphere(
                    dataset, ws.contains("name", "Carbon dioxide, fossil")
                ):
                    co2_amount = co2_flow["amount"]
                    co2_emitted = co2_amount * (1 - carbon_capture_rate)
                    co2_flow["amount"] = co2_emitted

                    # Add an input from the CCS dataset to the input dataset
                    ccs_exc = {
                        "uncertainty type": 0,
                        "loc": 0,
                        "amount": co2_amount - co2_emitted,
                        "type": "technosphere",
                        "production volume": 0,
                        "name": "carbon dioxide, captured at steel production plant, "
                        "with underground storage, post, 200 km",
                        "unit": "kilogram",
                        "location": dataset["location"],
                        "product": "carbon dioxide, captured and stored",
                    }
                    dataset["exchanges"].append(ccs_exc)

            if "log parameters" not in dataset:
                dataset["log parameters"] = {}

            dataset["log parameters"].update(
                {"carbon capture rate": carbon_capture_rate}
            )

        return datasets

    def write_log(self, dataset, status="created"):
        """
        Write log file.
        """

        logger.info(
            f"{status}|{self.model}|{self.scenario}|{self.year}|"
            f"{dataset['name']}|{dataset['location']}|"
            f"{dataset.get('log parameters', {}).get('carbon capture rate', '')}|"
            f"{dataset.get('log parameters', {}).get('thermal efficiency change', '')}|"
            f"{dataset.get('log parameters', {}).get('primary steel share', '')}|"
            f"{dataset.get('log parameters', {}).get('secondary steel share', '')}"
        )
