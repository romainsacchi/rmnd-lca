from .filesystem_constants import DATA_DIR
from .inventory_imports import DefaultInventory
from .transformation import BaseTransformation
from typing import List
from .data_collection import IAMDataCollection
from wurst import searching as ws

LCI_HEAT = DATA_DIR / "energy" / "lci-heat.xlsx"
FINAL_ENERGY_MAPPING = DATA_DIR / "energy" / "final_energy_mapping.xlsx"


class Energy(BaseTransformation):
    """
    Class that creates heating datasets based on IAM output data.

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
        cache: dict = None,
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
        )
        self.version = version

    def import_heating_inventories(
        self,
    ):
        """
        Import heating inventories into the database.

        :param database: The database to import into.
        :param system_model: The system model to import into.
        :param version: The version of the inventories to import.
        """

        inventory = DefaultInventory(
            database=self.database,
            version_in="3.9",
            version_out=self.version,
            path=LCI_HEAT,
            system_model=self.system_model,
            keep_uncertainty_data=False,
        )
        datasets = inventory.merge_inventory()
        self.database.extend(datasets)

        dataset_names = [d["name"] for d in datasets]
        self.regionalize_heating_datasets(dataset_names)

    def regionalize_heating_datasets(self, dataset_names: list):
        for ds in ws.get_many(
            self.database,
            ws.either(*[ws.equals("name", name) for name in dataset_names]),
            ws.equals("location", "RER"),
        ):
            new_datasets = self.fetch_proxies(
                ds["name"],
                ds["reference product"],
            )
            self.database.extend(new_datasets.values())
            self.add_to_index(new_datasets.values())
