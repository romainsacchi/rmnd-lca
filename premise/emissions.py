"""
Integrates projections regarding emissions of hot pollutants
from GAINS.
"""

from functools import lru_cache
from typing import Union

import numpy as np
import wurst
import xarray as xr
import yaml
from numpy import ndarray

from .filesystem_constants import DATA_DIR
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

logger = create_logger("emissions")

EI_POLLUTANTS = DATA_DIR / "GAINS_emission_factors" / "GAINS_ei_pollutants.yaml"
GAINS_SECTORS = DATA_DIR / "GAINS_emission_factors" / "GAINS_EU_sectors_mapping.yaml"


def fetch_mapping(filepath: str) -> dict:
    """Returns a dictionary from a YML file"""

    with open(filepath, "r", encoding="utf-8") as stream:
        mapping = yaml.safe_load(stream)
    return mapping


def _update_emissions(scenario, version, system_model, gains_scenario):
    emissions = Emissions(
        database=scenario["database"],
        year=scenario["year"],
        model=scenario["model"],
        pathway=scenario["pathway"],
        iam_data=scenario["iam data"],
        version=version,
        system_model=system_model,
        gains_scenario=gains_scenario,
    )

    emissions.update_emissions_in_database()
    scenario["database"] = emissions.database

    return scenario


class Emissions(BaseTransformation):
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
        gains_scenario: str,
    ):
        super().__init__(
            database,
            iam_data,
            model,
            pathway,
            year,
            version,
            system_model,
        )

        self.version = version
        self.gains_europe = self.prepare_data(iam_data.gains_data_EU)
        self.gains_global = self.prepare_data(iam_data.gains_data_IAM)
        self.ei_pollutants = fetch_mapping(EI_POLLUTANTS)
        self.gains_pollutant = {v: k for k, v in self.ei_pollutants.items()}
        self.gains_sectors = fetch_mapping(GAINS_SECTORS)
        self.gains_scenario = gains_scenario

        mapping = InventorySet(self.database)
        self.gains_map_europe: Dict[str, Set] = mapping.generate_gains_mapping()
        self.gains_map_global: Dict[str, Set] = mapping.generate_gains_mapping_IAM(
            mapping=self.gains_map_europe
        )
        self.rev_gains_map_europe, self.rev_gains_map_global = {}, {}

        for s in self.gains_map_europe:
            for t in self.gains_map_europe[s]:
                self.rev_gains_map_europe[t] = s

        for s in self.gains_map_global:
            for t in self.gains_map_global[s]:
                self.rev_gains_map_global[t] = s

    def prepare_data(self, data):

        def _(x):
            return xr.where((np.isnan(x)) | (x == 0), 1, x)

        if self.year in data.coords["year"].values:
            data = data.sel(year=[self.year]) / _(
                data.loc[
                    dict(
                        year=2020,
                    )
                ]
            )
        else:
            data = data.interp(year=[self.year]) / _(
                data.loc[
                    dict(
                        year=2020,
                    )
                ]
            )

        # replace 0 values with 1
        data = xr.where((np.isnan(data)) | (data == 0), 1, data)

        return data

    def update_emissions_in_database(self):
        # print("Integrating GAINS EU emission factors.")
        for ds in self.database:
            if (
                ds["name"] in self.rev_gains_map_europe
                and ds["location"] in self.gains_europe.coords["region"]
            ):
                gains_sector = self.rev_gains_map_europe[ds["name"]]
                self.update_pollutant_emissions(
                    ds,
                    gains_sector,
                    model="GAINS-EU",
                    regions=self.gains_europe.region.values,
                )
                self.write_log(ds, status="updated")

        # print("Integrating GAINS IAM emission factors.")
        for ds in self.database:
            if (
                ds["name"] in self.rev_gains_map_global
                and self.ecoinvent_to_iam_loc[ds["location"]]
                in self.gains_global.coords["region"]
            ):
                gains_sector = self.rev_gains_map_global[ds["name"]]
                self.update_pollutant_emissions(
                    ds,
                    gains_sector,
                    model="GAINS-IAM",
                    regions=self.gains_global.region.values,
                )

                self.write_log(ds, status="updated")

    def update_pollutant_emissions(
        self, dataset: dict, sector: str, model: str, regions: list
    ) -> dict:
        """
        Update pollutant emissions based on GAINS data.
        We apply a correction factor equal to the relative
        change in emissions compared to 2020

        :param dataset: dataset to adjust non-CO2 emission for
        :param sector: GAINS industrial sector to look up
        :return: Does not return anything. Modified in place.
        """

        # Update biosphere exchanges according to GAINS emission values
        for exc in ws.biosphere(
            dataset, ws.either(*[ws.equals("name", x) for x in self.ei_pollutants])
        ):
            gains_pollutant = self.ei_pollutants[exc["name"]]
            scaling_factor = self.find_gains_emissions_change(
                pollutant=gains_pollutant,
                location=(
                    dataset["location"]
                    if dataset["location"] in regions
                    else self.ecoinvent_to_iam_loc[dataset["location"]]
                ),
                sector=sector,
                model=model,
            )

            if scaling_factor != 1.0 and scaling_factor > 0.0:
                if f"{gains_pollutant} scaling factor" not in dataset.get(
                    "log parameters", {}
                ):
                    wurst.rescale_exchange(
                        exc, scaling_factor, remove_uncertainty=False
                    )

                    if "log parameters" not in dataset:
                        dataset["log parameters"] = {}

                    if "GAINS model" not in dataset["log parameters"]:
                        dataset["log parameters"]["GAINS model"] = model

                    if "GAINS sector" not in dataset["log parameters"]:
                        dataset["log parameters"]["GAINS sector"] = sector

                    dataset["log parameters"].update(
                        {f"{gains_pollutant} scaling factor": scaling_factor}
                    )

        return dataset

    @lru_cache
    def find_gains_emissions_change(
        self, pollutant: str, location: str, sector: str, model: str
    ) -> Union[ndarray, float]:
        """
        Return the relative change in emissions compared to 2020
        for a given pollutant, location and sector.
        :param pollutant: name of pollutant
        :param sector: name of technology/sector
        :param location: location of emitting dataset
        :model: GAINS model
        :return: a scaling factor
        """

        data = self.gains_europe if model == "GAINS-EU" else self.gains_global

        key_exists = all(
            k in data.coords[dim].values
            for k, dim in zip(
                [location, pollutant, sector, self.year],
                ["region", "pollutant", "sector", "year"],
            )
        )

        if key_exists:
            scaling_factor = data.loc[
                dict(
                    region=location,
                    pollutant=pollutant,
                    sector=sector,
                )
            ]

            scaling_factor = np.clip(
                scaling_factor,
                1 if self.year < 2020 else 0,
                1 if self.year > 2020 else 1e6,
            )

            if np.isnan(scaling_factor) or scaling_factor == 0.0:
                scaling_factor = 1.0

            return float(scaling_factor)
        return 1.0

    def write_log(self, dataset, status="created"):
        """
        Write log file.
        """

        if "log parameters" in dataset:
            if "GAINS model" in dataset["log parameters"]:
                logger.info(
                    f"{status}|{self.model}|{self.scenario}|{self.year}|"
                    f"{dataset['name']}|{dataset['location']}|"
                    f"{dataset.get('log parameters', {}).get('GAINS model', '')}|"
                    f"{dataset.get('log parameters', {}).get('GAINS sector', '')}|"
                    f"{dataset.get('log parameters', {}).get('CO scaling factor', '')}|"
                    f"{dataset.get('log parameters', {}).get('CH4 scaling factor', '')}|"
                    f"{dataset.get('log parameters', {}).get('N2O scaling factor', '')}|"
                    f"{dataset.get('log parameters', {}).get('NH3 scaling factor', '')}|"
                    f"{dataset.get('log parameters', {}).get('NOx scaling factor', '')}|"
                    f"{dataset.get('log parameters', {}).get('PM1 scaling factor', '')}|"
                    f"{dataset.get('log parameters', {}).get('PM10 scaling factor', '')}|"
                    f"{dataset.get('log parameters', {}).get('PM25 scaling factor', '')}|"
                    f"{dataset.get('log parameters', {}).get('SO2 scaling factor', '')}|"
                    f"{dataset.get('log parameters', {}).get('VOC scaling factor', '')}"
                )
