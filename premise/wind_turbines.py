"""
module to create a variant of the onshore and offshore wind
turbine datasets to represent the direct-drive technology.

"""

import copy
import uuid

from .activity_maps import InventorySet
from .logger import create_logger
from .transformation import BaseTransformation, IAMDataCollection, List, np, ws

logger = create_logger("wind_turbine")


def _update_wind_turbines(scenario, version, system_model):
    wind_turbine = WindTurbine(
        database=scenario["database"],
        iam_data=scenario["iam data"],
        model=scenario["model"],
        pathway=scenario["pathway"],
        year=scenario["year"],
        version=version,
        system_model=system_model,
        cache=scenario.get("cache"),
        index=scenario.get("index"),
    )

    wind_turbine.create_direct_drive_turbines()

    scenario["database"] = wind_turbine.database
    scenario["index"] = wind_turbine.index
    scenario["cache"] = wind_turbine.cache

    return scenario


class WindTurbine(BaseTransformation):
    """
    Class that create additional wind turbine datasets
    to represent the direct-drive technology.

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
        cache: dict = None,
        index: dict = None,
    ) -> None:
        super().__init__(
            database,
            iam_data,
            model,
            pathway,
            year,
            version,
            system_model,
            cache,
            index,
        )
        self.system_model = system_model
        mapping = InventorySet(database=database, version=version, model=model)
        self.powerplant_map = mapping.generate_powerplant_map()

    def create_direct_drive_turbines(self):
        """
        Create direct-drive wind turbine datasets.
        """

        wind_turbine_technologies = ["Wind Onshore", "Wind Offshore"]

        for technology in wind_turbine_technologies:
            for dataset in ws.get_many(
                self.database,
                ws.either(
                    *[
                        ws.equals("name", tech)
                        for tech in self.powerplant_map[technology]
                    ]
                ),
            ):
                dataset_copy = self.create_dataset_copy(dataset, "direct drive")
                self.database.append(dataset_copy)

                construction_datasets = self.create_construction_datasets()
                self.database.extend(construction_datasets)

                market_datasets = self.create_market_datasets(
                    dataset_copy, construction_datasets
                )
                self.database.extend(market_datasets)

                self.update_production_dataset_links(dataset_copy, market_datasets)

    def create_dataset_copy(self, dataset, suffix):
        dataset_copy = copy.deepcopy(dataset)
        dataset_copy["name"] += f", {suffix}"
        dataset_copy["code"] = str(uuid.uuid4().hex)
        dataset_copy["comment"] = (
            f"This dataset represents the {suffix} technology "
            "variant of the ecoinvent dataset."
        )

        for exc in ws.production(dataset_copy):
            exc["name"] = dataset_copy["name"]
            if "input" in exc:
                del exc["input"]

        self.write_log(dataset_copy)
        return dataset_copy

    def create_construction_datasets(self):
        construction_names = [
            "wind power plant construction, 800kW, moving parts",
            "wind power plant construction, 2MW, offshore, moving parts",
            "wind turbine construction, 4.5MW, onshore",
            "wind turbine construction, 2MW, onshore",
        ]

        construction_datasets = []
        for name in construction_names:
            construction_dataset = ws.get_one(self.database, ws.equals("name", name))
            construction_dataset_copy = self.create_dataset_copy(
                construction_dataset, "direct drive"
            )
            construction_datasets.append(construction_dataset_copy)

        return construction_datasets

    def create_market_datasets(self, parent_dataset, construction_datasets):
        market_names = [
            "market for wind power plant, 800kW, moving parts",
            "market for wind power plant, 2MW, offshore, moving parts",
            "market for wind turbine, 4.5MW, onshore",
            "market for wind turbine, 2MW, onshore",
        ]

        product_names = [
            "wind power plant, 800kW, moving parts",
            "wind power plant, 2MW, offshore, moving parts",
            "wind turbine, 4.5MW, onshore",
            "wind turbine, 2MW, onshore",
        ]

        market_datasets = []
        for market_name, product_name, construction_dataset in zip(
            market_names, product_names, construction_datasets
        ):
            market_dataset = ws.get_one(self.database, ws.equals("name", market_name))
            market_dataset_copy = self.create_dataset_copy(
                market_dataset, "direct drive"
            )

            for exc in ws.technosphere(market_dataset_copy):
                if exc["product"] == product_name:
                    exc["name"] = construction_dataset["name"]
                    exc["product"] = product_name

            self.write_log(market_dataset_copy)
            market_datasets.append(market_dataset_copy)

        return market_datasets

    @staticmethod
    def update_production_dataset_links(production_dataset, market_datasets):
        for market_dataset in market_datasets:
            market_product = [exc for exc in ws.production(market_dataset)][0][
                "product"
            ]
            market_name = market_dataset["name"]
            for exc in ws.technosphere(production_dataset):
                if exc["product"] == market_product:
                    exc["name"] = market_name
                    exc["product"] = market_product

    def write_log(self, dataset, status="created"):
        """
        Write log file.
        """

        logger.info(
            f"{status}|{self.model}|{self.scenario}|{self.year}|"
            f"{dataset['name']}|{dataset['location']}"
        )
