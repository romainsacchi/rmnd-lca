"""
Integrates projections regarding direct air capture and storage.
"""

import copy

import yaml

from .filesystem_constants import DATA_DIR
from .logger import create_logger
from .transformation import (
    BaseTransformation,
    IAMDataCollection,
    InventorySet,
    List,
    uuid,
    ws,
)
from .utils import rescale_exchanges

logger = create_logger("dac")

HEAT_SOURCES = DATA_DIR / "fuels" / "heat_sources_map.yml"


def fetch_mapping(filepath: str) -> dict:
    """Returns a dictionary from a YML file"""

    with open(filepath, "r", encoding="utf-8") as stream:
        mapping = yaml.safe_load(stream)
    return mapping


def _update_dac(scenario, version, system_model):
    dac = DirectAirCapture(
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

    if scenario["iam data"].dac_markets is not None:
        dac.generate_dac_activities()
        dac.relink_datasets()
        scenario["database"] = dac.database
        scenario["cache"] = dac.cache
        scenario["index"] = dac.index
    else:
        print("No DAC information found in IAM data. Skipping.")

    return scenario


class DirectAirCapture(BaseTransformation):
    """
    Class that modifies DAC and DACCS inventories and markets
    in ecoinvent based on IAM output data.
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
    ):
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
        self.database = database
        self.iam_data = iam_data
        self.model = model
        self.pathway = pathway
        self.year = year
        self.version = version
        self.system_model = system_model
        mapping = InventorySet(self.database)
        self.dac_plants = mapping.generate_daccs_map()
        self.carbon_storage = mapping.generate_carbon_storage_map()

    def generate_dac_activities(self) -> None:
        """
        Generates regional variants of the direct air capture process with varying heat sources.

        This function fetches the original datasets for the direct air capture process and creates regional variants
        with different heat sources. The function loops through the heat sources defined in the `HEAT_SOURCES` mapping,
        modifies the original datasets to include the heat source, and adds the modified datasets to the database.

        """
        # print("Generate region-specific direct air capture processes.")

        # get original dataset
        for ds_list in self.carbon_storage.values():
            for ds_name in ds_list:
                new_ds = self.fetch_proxies(
                    name=ds_name,
                    ref_prod="carbon dioxide, stored",
                )

                for k, dataset in new_ds.items():
                    # Add created dataset to cache
                    self.add_new_entry_to_cache(
                        location=dataset["location"],
                        exchange=dataset,
                        allocated=[dataset],
                        shares=[
                            1.0,
                        ],
                    )
                    self.write_log(dataset)
                    # add it to list of created datasets
                    self.add_to_index(dataset)

                self.database.extend(new_ds.values())

        # define heat sources
        heat_map_ds = fetch_mapping(HEAT_SOURCES)

        # get original dataset
        for technology, ds_list in self.dac_plants.items():
            for ds_name in ds_list:
                original_ds = self.fetch_proxies(
                    name=ds_name,
                    ref_prod="carbon dioxide",
                    relink=False,
                    delete_original_dataset=False,
                    empty_original_activity=False,
                )

                # loop through heat sources
                for heat_type, activities in heat_map_ds.items():
                    # with consequential modeling, waste heat is not available
                    if (
                        self.system_model == "consequential"
                        and heat_type == "waste heat"
                    ):
                        continue

                    # with liquid solvent-based DAC, we cannot use waste heat
                    # because the operational temperature required is 900C
                    if technology in ["dac_solvent", "daccs_solvent"]:
                        if heat_type == "waste heat":
                            continue

                    new_ds = copy.deepcopy(original_ds)
                    for k, dataset in new_ds.items():
                        dataset["name"] += f", with {heat_type}, and grid electricity"
                        dataset["code"] = str(uuid.uuid4().hex)
                        dataset["comment"] += activities["description"]

                        for exc in ws.production(dataset):
                            exc["name"] = dataset["name"]
                            if "input" in exc:
                                del exc["input"]

                        for exc in ws.technosphere(dataset):
                            if "heat" in exc["name"]:
                                exc["name"] = activities["name"]
                                exc["product"] = activities["reference product"]
                                exc["location"] = "RoW"

                                if "input" in exc:
                                    del exc["input"]

                                if heat_type == "heat pump heat":
                                    exc["unit"] = "kilowatt hour"
                                    exc["location"] = "GLO"
                                    exc["amount"] *= 1 / (2.9 * 3.6)

                        new_ds[k] = self.relink_technosphere_exchanges(
                            dataset,
                        )

                    # adjust efficiency, if needed
                    new_ds = self.adjust_dac_efficiency(new_ds)

                    self.database.extend(new_ds.values())

                    # add to log
                    for dataset in list(new_ds.values()):
                        self.write_log(dataset)
                        # add it to list of created datasets
                        self.add_to_index(dataset)

    def adjust_dac_efficiency(self, datasets):
        """
        Fetch the cumulated deployment of DAC from IAM file.
        Apply a learning rate -- see Qiu et al., 2022.
        """

        for region, dataset in datasets.items():
            if self.iam_data.dac_electricity_efficiencies is not None:
                if (
                    region
                    in self.iam_data.dac_electricity_efficiencies.coords[
                        "region"
                    ].values
                ):
                    if (
                        self.year
                        in self.iam_data.dac_electricity_efficiencies.coords[
                            "year"
                        ].values
                    ):
                        scaling_factor = float(
                            1
                            / self.iam_data.dac_electricity_efficiencies.sel(
                                region=region, year=self.year
                            ).values.item()
                        )
                    else:
                        scaling_factor = float(
                            1
                            / self.iam_data.dac_electricity_efficiencies.sel(
                                region=region
                            )
                            .interp(year=self.year)
                            .values
                        )

                    # bound the scaling factor to 1.5 and 0.5
                    scaling_factor = max(0.5, min(1.5, scaling_factor))

                    if scaling_factor != 1:
                        rescale_exchanges(
                            dataset,
                            scaling_factor,
                            technosphere_filters=[ws.equals("unit", "kilowatt hour")],
                        )

                        # add in comments the scaling factor applied
                        dataset["comment"] += (
                            f" The electrical efficiency of the system has been "
                            f"adjusted to match the efficiency of the "
                            f"average DAC plant in {self.year}."
                        )

                        if "log parameters" not in dataset:
                            dataset["log parameters"] = {}

                        dataset["log parameters"].update(
                            {
                                "electricity scaling factor": scaling_factor,
                            }
                        )

            if self.iam_data.dac_heat_efficiencies is not None:
                if (
                    region
                    in self.iam_data.dac_heat_efficiencies.coords["region"].values
                ):
                    if (
                        self.year
                        in self.iam_data.dac_heat_efficiencies.coords["year"].values
                    ):
                        scaling_factor = float(
                            1
                            / self.iam_data.dac_heat_efficiencies.sel(
                                region=region, year=self.year
                            ).values.item()
                        )
                    else:
                        scaling_factor = float(
                            1
                            / self.iam_data.dac_heat_efficiencies.sel(region=region)
                            .interp(year=self.year)
                            .values
                        )

                    # bound the scaling factor to 1.5 and 0.5
                    scaling_factor = max(0.5, min(1.5, scaling_factor))

                    if scaling_factor != 1:

                        rescale_exchanges(
                            dataset,
                            scaling_factor,
                            technosphere_filters=[ws.equals("unit", "megajoule")],
                        )

                        # add in comments the scaling factor applied
                        dataset["comment"] += (
                            f" The thermal efficiency of the system has been "
                            f"adjusted to match the efficiency of the "
                            f"average DAC plant in {self.year}."
                        )

                        if "log parameters" not in dataset:
                            dataset["log parameters"] = {}

                        dataset["log parameters"].update(
                            {
                                "heat scaling factor": scaling_factor,
                            }
                        )

        return datasets

    def write_log(self, dataset, status="created"):
        """
        Write log file.
        """
        logger.info(
            f"{status}|{self.model}|{self.scenario}|{self.year}|"
            f"{dataset['name']}|{dataset['location']}|"
            f"{dataset.get('log parameters', {}).get('electricity scaling factor', '')}|"
            f"{dataset.get('log parameters', {}).get('heat scaling factor', '')}"
        )
