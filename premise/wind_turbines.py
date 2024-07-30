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


def relink(
    dataset,
) -> dict:
    """
    Relink technosphere exchanges to the new datasets.
    """

    for exc in ws.technosphere(
        dataset,
        ws.equals("unit", "unit"),
        ws.exclude(ws.contains("name", "connection")),
    ):
        exc["name"] += ", direct drive"

    return dataset


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

        datasets_terms = [
            "electricity production, wind, <1MW turbine, onshore",
            "electricity production, wind, >3MW turbine, onshore",
            "electricity production, wind, 1-3MW turbine, onshore",
            "electricity production, wind, 1-3MW turbine, offshore",
            "wind power plant construction",
            "wind turbine construction",
            "market for wind power plant",
            "market for wind turbine",
        ]

        new_datasets, processed = [], []
        for dataset in ws.get_many(
            self.database,
            ws.either(*[ws.contains("name", tech) for tech in datasets_terms]),
            ws.exclude(ws.contains("name", "direct drive")),
        ):
            if dataset["code"] not in processed:
                dataset_copy = self.create_dataset_copy(dataset, "direct drive")
                dataset_copy = relink(dataset_copy)
                new_datasets.append(dataset_copy)
                processed.append(dataset["code"])

        self.database.extend(new_datasets)

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

    def write_log(self, dataset, status="created"):
        """
        Write log file.
        """

        logger.info(
            f"{status}|{self.model}|{self.scenario}|{self.year}|"
            f"{dataset['name']}|{dataset['location']}"
        )
