"""
Integrates projections regarding direct air capture and storage.
"""

import copy
import logging.config
from pathlib import Path

import wurst
import yaml

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

logger = logging.getLogger("dac")

import numpy as np

from .transformation import (
    BaseTransformation,
    IAMDataCollection,
    InventorySet,
    List,
    uuid,
    ws,
)

HEAT_SOURCES = DATA_DIR / "fuels" / "heat_sources_map.yml"


def fetch_mapping(filepath: str) -> dict:
    """Returns a dictionary from a YML file"""

    with open(filepath, "r", encoding="utf-8") as stream:
        mapping = yaml.safe_load(stream)
    return mapping


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
        print("Generate region-specific direct air capture processes.")

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
                    delete_original_dataset=True,
                )

                # loop through heat sources
                for heat_type, activities in heat_map_ds.items():
                    # with consequential modeling, waste heat is not available
                    if (
                        self.system_model == "consequential"
                        and heat_type == "waste heat"
                    ):
                        continue

                    # with solvent-based DAC, we cannot use waste heat
                    # because the operational temperature required is 900C
                    if technology in ["solvent_dac", "solvent_daccs"]:
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
                                    exc["amount"] *= 1 / (2.9 * 3.6)

                        new_ds[k] = self.relink_technosphere_exchanges(
                            dataset,
                        )

                    # adjust efficiency, if needed
                    new_ds = self.adjust_dac_efficiency(new_ds, technology)

                    self.database.extend(new_ds.values())

                    # add to log
                    for dataset in list(new_ds.values()):
                        self.write_log(dataset)
                        # add it to list of created datasets
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

    def adjust_dac_efficiency(self, datasets, technology):
        """
        Fetch the cumulated deployment of DAC from IAM file.
        Apply a learning rate -- see Qiu et al., 2022.
        """

        # learning rates for operation-related expenditures
        # (thermal and electrical energy)
        learning_rates_operation = {
            "dac_solvent": 0.025,
            "dac_sorbent": 0.025,
            "daccs_solvent": 0.025,
            "daccs_sorbent": 0.025,
        }

        # learning rates for
        # infrastructure-related expenditures
        learning_rates_infra = {
            "dac_solvent": 0.1,
            "dac_sorbent": 0.15,
            "daccs_solvent": 0.1,
            "daccs_sorbent": 0.15,
        }

        theoretical_min_operation = {
            "dac_solvent": 0.95,
            "dac_sorbent": 0.95,
            "daccs_solvent": 0.95,
            "daccs_sorbent": 0.95,
        }

        theoretical_min_infra = {
            "dac_solvent": 0.44,
            "dac_sorbent": 0.18,
            "daccs_solvent": 0.44,
            "daccs_sorbent": 0.18,
        }

        for region, dataset in datasets.items():
            # fetch cumulated deployment of DAC from IAM file
            if "dac_solvent" in self.iam_data.production_volumes.variables.values:
                cumulated_deployment = (
                    np.clip(
                        self.iam_data.production_volumes.sel(
                            variables="dac_solvent",
                        )
                        .interp(year=self.year)
                        .sum(dim="region")
                        .values.item()
                        * -1,
                        1e-3,
                        None,
                    )
                    / 2
                )  # divide by 2,
                # as we assume sorbent and solvent are deployed equally

                initial_deployment = (
                    np.clip(
                        self.iam_data.production_volumes.sel(
                            variables="dac_solvent", year=2020
                        )
                        .sum(dim="region")
                        .values.item()
                        * -1,
                        1e-3,
                        None,
                    )
                    / 2
                )  # divide by 2,
                # as we assume sorbent and solvent are deployed equally
            else:
                cumulated_deployment = 1
                initial_deployment = 1

            # the learning rate is applied per doubling
            # of the cumulative deployment
            # relative to 2020

            scaling_factor_operation = (
                1 - theoretical_min_operation[technology]
            ) * np.power(
                (1 - learning_rates_operation[technology]),
                np.log2(cumulated_deployment / initial_deployment),
            ) + theoretical_min_operation[
                technology
            ]

            scaling_factor_infra = (1 - theoretical_min_infra[technology]) * np.power(
                (1 - learning_rates_infra[technology]),
                np.log2(cumulated_deployment / initial_deployment),
            ) + theoretical_min_infra[technology]

            current_energy_inputs = sum(
                e["amount"] for e in dataset["exchanges"] if e["unit"] == "megajoule"
            )
            current_energy_inputs += sum(
                e["amount"] * 3.6
                for e in dataset["exchanges"]
                if e["unit"] == "kilowatt hour"
            )

            if "log parameters" not in dataset:
                dataset["log parameters"] = {}

            dataset["log parameters"].update(
                {
                    "initial energy input per kg CO2": current_energy_inputs,
                }
            )

            if scaling_factor_operation != 1:
                # Scale down the energy exchanges using the scaling factor
                wurst.change_exchanges_by_constant_factor(
                    dataset,
                    scaling_factor_operation,
                    technosphere_filters=[
                        ws.either(
                            *[ws.contains("name", x) for x in ["heat", "electricity"]]
                        )
                    ],
                    biosphere_filters=[ws.exclude(ws.contains("type", "biosphere"))],
                )

            new_energy_inputs = sum(
                e["amount"] for e in dataset["exchanges"] if e["unit"] == "megajoule"
            )
            new_energy_inputs += sum(
                e["amount"] * 3.6
                for e in dataset["exchanges"]
                if e["unit"] == "kilowatt hour"
            )

            dataset["log parameters"].update(
                {
                    "new energy input per kg CO2": new_energy_inputs,
                }
            )

            # add in comments the scaling factor applied
            dataset["comment"] += (
                f" Operation-related expenditures have been "
                f"reduced by: {int((1 - scaling_factor_operation) * 100)}%."
            )

            dataset["log parameters"].update(
                {
                    "scaling factor operation": scaling_factor_operation,
                    "scaling factor infrastructure": scaling_factor_infra,
                }
            )

            if scaling_factor_infra != 1:
                # Scale down the infra and material exchanges using the scaling factor
                wurst.change_exchanges_by_constant_factor(
                    dataset,
                    scaling_factor_infra,
                    technosphere_filters=[
                        ws.exclude(
                            ws.either(
                                *[
                                    ws.contains("name", x)
                                    for x in ["heat", "electricity", "storage"]
                                ]
                            )
                        )
                    ],
                    biosphere_filters=[ws.exclude(ws.contains("type", "biosphere"))],
                )

            # add in comments the scaling factor applied
            dataset["comment"] += (
                f" Infrastructure-related expenditures have been "
                f"reduced by: {int((1 - scaling_factor_infra) * 100)}%."
            )

        return datasets

    def write_log(self, dataset, status="created"):
        """
        Write log file.
        """
        logger.info(
            f"{status}|{self.model}|{self.scenario}|{self.year}|"
            f"{dataset['name']}|{dataset['location']}|"
            f"{dataset.get('log parameters', {}).get('scaling factor operation', '')}|"
            f"{dataset.get('log parameters', {}).get('scaling factor infrastructure', '')}|"
            f"{dataset.get('log parameters', {}).get('initial energy input per kg CO2', '')}|"
            f"{dataset.get('log parameters', {}).get('new energy input per kg CO2', '')}"
        )
