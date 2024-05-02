"""
transport.py contains the class Transport, which takes care of importing inventories
for a number of different vehicle types, and create fleet average vehicles based on
IAM data, and integrate them into the database.
"""

import re
import uuid
from typing import Any, Dict, List, Union

import numpy as np
import xarray as xr
import yaml
from wurst import searching as ws
from wurst import transformations as wt

from .filesystem_constants import DATA_DIR, IAM_OUTPUT_DIR, INVENTORY_DIR
from .inventory_imports import VariousVehicles
from .transformation import BaseTransformation, IAMDataCollection
from .utils import HiddenPrints, eidb_label
from .validation import CarValidation, TruckValidation

FILEPATH_FLEET_COMP = IAM_OUTPUT_DIR / "fleet_files" / "fleet_all_vehicles.csv"
FILEPATH_IMAGE_TRUCKS_FLEET_COMP = (
    IAM_OUTPUT_DIR / "fleet_files" / "image_fleet_trucks.csv"
)
FILEPATH_TRUCK_LOAD_FACTORS = DATA_DIR / "transport" / "avg_load_factors.yaml"
FILEPATH_VEHICLES_MAP = DATA_DIR / "transport" / "vehicles_map.yaml"


def _update_vehicles(scenario, vehicle_type, version, system_model):

    has_fleet = False
    if vehicle_type == "car":
        if hasattr(scenario["iam data"], "trsp_cars"):
            has_fleet = True
    elif vehicle_type == "truck":
        if hasattr(scenario["iam data"], "trsp_trucks"):
            has_fleet = True
    elif vehicle_type == "bus":
        if hasattr(scenario["iam data"], "trsp_buses"):
            has_fleet = True
    elif vehicle_type == "two wheeler":
        if hasattr(scenario["iam data"], "trsp_two_wheelers"):
            has_fleet = True
    else:
        raise ValueError("Unknown vehicle type.")

    trspt = Transport(
        database=scenario["database"],
        year=scenario["year"],
        model=scenario["model"],
        pathway=scenario["pathway"],
        iam_data=scenario["iam data"],
        version=version,
        system_model=system_model,
        vehicle_type=vehicle_type,
        relink=False,
        has_fleet=has_fleet,
        index=scenario.get("index"),
    )

    scenario["database"] = trspt.database
    scenario["cache"] = trspt.cache
    scenario["index"] = trspt.index

    if vehicle_type == "car":
        validate = CarValidation(
            model=scenario["model"],
            scenario=scenario["pathway"],
            year=scenario["year"],
            regions=scenario["iam data"].regions,
            database=trspt.database,
            iam_data=scenario["iam data"],
        )

        validate.run_car_checks()

    if vehicle_type == "truck":
        validate = TruckValidation(
            model=scenario["model"],
            scenario=scenario["pathway"],
            year=scenario["year"],
            regions=scenario["iam data"].regions,
            database=trspt.database,
            iam_data=scenario["iam data"],
        )

        validate.run_truck_checks()

    return scenario


def get_average_truck_load_factors() -> Dict[str, Dict[str, Dict[str, float]]]:
    """
    Load average load factors for trucks
    to convert transport demand in vkm into tkm.
    :return: dictionary with load factors per truck size class
    """
    with open(FILEPATH_TRUCK_LOAD_FACTORS, "r", encoding="utf-8") as stream:
        out = yaml.safe_load(stream)
        return out


def get_vehicles_mapping() -> Dict[str, dict]:
    """
    Return a dictionary that contains mapping
    between `carculator` terminology and `ecoinvent` terminology
    regarding size classes, powertrain types, etc.
    :return: dictionary to map terminology between carculator and ecoinvent
    """
    with open(FILEPATH_VEHICLES_MAP, "r", encoding="utf-8") as stream:
        out = yaml.safe_load(stream)
        return out


def normalize_exchange_amounts(list_act: List[dict]) -> List[dict]:
    """
    In vehicle market datasets, we need to ensure that the total contribution
    of single vehicle types equal 1.

    :param list_act: list of transport market activities
    :return: same list, with activity exchanges normalized to 1

    """

    for act in list_act:
        total = 0
        for exc in act["exchanges"]:
            if exc["type"] == "technosphere":
                total += exc["amount"]

        for exc in act["exchanges"]:
            if exc["type"] == "technosphere":
                exc["amount"] /= total

    return list_act


def create_fleet_vehicles(
    datasets: List[dict],
    vehicle_type: str,
    year: int,
    model: str,
    scenario: str,
    version: str,
    system_model: str,
    regions: List[str],
    arr: xr.DataArray,
) -> List[dict[str, Union[Union[str, float], Any]]]:
    """
    Create datasets for fleet average vehicles based on IAM fleet data.

    :param datasets: vehicle datasets of all size, powertrain and construction years.
    :param vehicle_type: "car", "truck"
    :param year: year for the fleet average vehicle
    :param model: IAM model
    :param scenario: IAM scenario
    :param regions: IAM regions
    :return: list of fleet average vehicle datasets
    """

    vehicles_map = get_vehicles_mapping()

    list_act = []

    # average load factors for trucks
    # to convert vkm to tkm
    avg_load = get_average_truck_load_factors()

    # add missing IMAGE regions
    d_missing_regions = {
        "BRA": "RSAM",
        "CEU": "WEU",
        "CAN": "OCE",
        "KOR": "SEAS",
        "SAF": "WAF",
        "RUS": "UKR",
        "INDO": "SEAS",
        "ME": "TUR",
        "RSAF": "WAF",
        "EAF": "WAF",
        "MEX": "RSAM",
        "NAF": "WAF",
        "RCAM": "RSAM",
        "RSAS": "SEAS",
        "STAN": "TUR",
    }

    for region in regions:
        fleet_region = d_missing_regions.get(region, region)

        sizes = [
            s
            for s in vehicles_map[vehicle_type]["sizes"]
            if s in arr.coords["size"].values
        ]

        if year in arr.coords["year"].values:
            region_size_fleet = arr.sel(
                region=fleet_region, size=sizes, year=int(np.clip(year, 2005, 2050))
            )
        else:
            region_size_fleet = arr.sel(region=fleet_region, size=sizes).interp(
                year=int(np.clip(year, 2005, 2050))
            )

        total_km = region_size_fleet.sum()

        if total_km > 0:
            name = (
                f"{vehicles_map[vehicle_type]['name']}, unspecified, long haul"
                if vehicle_type == "truck"
                else f"{vehicles_map[vehicle_type]['name']}, unspecified"
            )
            act = {
                "name": name,
                "reference product": vehicles_map[vehicle_type]["name"],
                "unit": vehicles_map[vehicle_type]["unit"],
                "location": region,
                "exchanges": [
                    {
                        "name": name,
                        "product": vehicles_map[vehicle_type]["name"],
                        "unit": vehicles_map[vehicle_type]["unit"],
                        "location": region,
                        "type": "production",
                        "amount": 1,
                    }
                ],
                "code": str(uuid.uuid4().hex),
                "database": "premise",
                "comment": f"Fleet-average vehicle for the year {year}, "
                f"for the region {region}.",
            }

            for size in sizes:
                for pwt in arr.coords["powertrain"].values:
                    indiv_km = region_size_fleet.sel(
                        size=size,
                        powertrain=pwt,
                    )
                    if indiv_km > 0:
                        indiv_share = (indiv_km / total_km).values.item(0)

                        load = 1
                        if vehicle_type == "truck":
                            load = avg_load[vehicle_type]["long haul"][size]

                        filters = [
                            ws.contains("name", size),
                            ws.contains("name", pwt),
                            ws.equals("location", region),
                        ]

                        if vehicle_type != "truck":
                            filters.append(ws.contains("name", vehicle_type))

                        if pwt in ["diesel", "gasoline", "compressed gas"]:
                            filters.extend(
                                [
                                    ws.exclude(ws.contains("name", "plugin")),
                                    ws.exclude(ws.contains("name", "hybrid")),
                                    ws.exclude(ws.contains("name", "ab")),
                                ]
                            )

                        if pwt in ["diesel hybrid", "gasoline hybrid"]:
                            filters.extend(
                                [
                                    ws.exclude(ws.contains("name", "plugin")),
                                    ws.exclude(ws.contains("name", "ab")),
                                ]
                            )

                        if pwt in ["plugin diesel hybrid", "plugin gasoline hybrid"]:
                            filters.extend(
                                [
                                    ws.exclude(ws.contains("name", "ab")),
                                ]
                            )

                        if pwt not in [
                            "battery electric",
                            "fuel cell electric",
                            "battery electric - overnight charging",
                        ]:
                            if vehicle_type == "car":
                                filters.append(ws.contains("name", "EURO-6"))
                            if vehicle_type in ["bus", "truck"]:
                                filters.append(ws.contains("name", "EURO-VI"))

                        if size in ["Medium", "Large"]:
                            filters.append(ws.exclude(ws.contains("name", "SUV")))

                        try:
                            vehicle = ws.get_one(datasets, *filters)

                        except ws.MultipleResults:
                            print(
                                "Multiple results for", size, vehicle_type, pwt, region
                            )
                            vehicles = ws.get_many(datasets, *filters)
                            for vehicle in vehicles:
                                print(vehicle["name"], vehicle["location"])
                            raise
                        except ws.NoResults:
                            print("No results for", size, vehicle_type, pwt, region)
                            vehicle = {
                                "name": "unknown",
                                "reference product": "unknown",
                                "unit": "unknown",
                            }

                        act["exchanges"].append(
                            {
                                "name": vehicle["name"],
                                "product": vehicle["reference product"],
                                "unit": vehicle["unit"],
                                "location": region,
                                "type": "technosphere",
                                "amount": indiv_share * load,
                            }
                        )

            if len(act["exchanges"]) > 1:
                list_act.append(act)

            # also create size-specific fleet vehicles
            if vehicle_type == "truck":
                for size in sizes:
                    total_size_km = region_size_fleet.sel(size=size).sum()

                    if total_size_km > 0:
                        name = (
                            f"{vehicles_map[vehicle_type]['name']}, {size} gross weight, "
                            f"unspecified powertrain, long haul"
                        )
                        act = {
                            "name": name,
                            "reference product": vehicles_map[vehicle_type]["name"],
                            "unit": vehicles_map[vehicle_type]["unit"],
                            "location": region,
                            "exchanges": [
                                {
                                    "name": name,
                                    "product": vehicles_map[vehicle_type]["name"],
                                    "unit": vehicles_map[vehicle_type]["unit"],
                                    "location": region,
                                    "type": "production",
                                    "amount": 1,
                                }
                            ],
                            "code": str(uuid.uuid4().hex),
                            "database": eidb_label(
                                {"model": model, "pathway": scenario, "year": year},
                                version,
                                system_model,
                            ),
                            "comment": f"Fleet-average vehicle for the year {year}, for the region {region}.",
                        }

                        for pwt in region_size_fleet.coords["powertrain"].values:
                            indiv_km = region_size_fleet.sel(
                                size=size,
                                powertrain=pwt,
                            )
                            if indiv_km > 0:
                                indiv_share = (indiv_km / total_size_km).values.item(0)
                                load = avg_load[vehicle_type]["long haul"][size]

                                filters = [
                                    ws.contains("name", size),
                                    ws.contains("name", pwt),
                                    ws.equals("location", region),
                                ]

                                if vehicle_type != "truck":
                                    filters.append(ws.contains("name", vehicle_type))

                                if pwt in ["diesel", "gasoline", "compressed gas"]:
                                    filters.extend(
                                        [
                                            ws.exclude(ws.contains("name", "plugin")),
                                            ws.exclude(ws.contains("name", "hybrid")),
                                            ws.exclude(ws.contains("name", "ab")),
                                        ]
                                    )

                                if pwt in ["diesel hybrid", "gasoline hybrid"]:
                                    filters.extend(
                                        [
                                            ws.exclude(ws.contains("name", "plugin")),
                                            ws.exclude(ws.contains("name", "ab")),
                                        ]
                                    )

                                if pwt in [
                                    "plugin diesel hybrid",
                                    "plugin gasoline hybrid",
                                ]:
                                    filters.extend(
                                        [
                                            ws.exclude(ws.contains("name", "ab")),
                                        ]
                                    )

                                if pwt not in [
                                    "battery electric",
                                    "fuel cell electric",
                                    "battery electric - overnight charging",
                                ]:
                                    if vehicle_type in ["bus", "truck"]:
                                        filters.append(ws.contains("name", "EURO-VI"))

                                try:
                                    vehicle = ws.get_one(datasets, *filters)

                                except ws.MultipleResults:
                                    print(
                                        "Multiple results for",
                                        size,
                                        vehicle_type,
                                        pwt,
                                        region,
                                    )
                                    vehicles = ws.get_many(datasets, *filters)
                                    for vehicle in vehicles:
                                        print(vehicle["name"], vehicle["location"])
                                    raise

                                except ws.NoResults:
                                    print(
                                        "No results for SIZE",
                                        size,
                                        vehicle_type,
                                        pwt,
                                        region,
                                    )
                                    vehicle = {
                                        "name": "unknown",
                                        "reference product": "unknown",
                                        "unit": "unknown",
                                    }

                                act["exchanges"].append(
                                    {
                                        "name": vehicle["name"],
                                        "product": vehicle["reference product"],
                                        "unit": vehicle["unit"],
                                        "location": region,
                                        "type": "technosphere",
                                        "amount": indiv_share * load,
                                    }
                                )

                        if len(act["exchanges"]) > 1:
                            list_act.append(act)

    return normalize_exchange_amounts(list_act)


class Transport(BaseTransformation):
    """
    Class that modifies transport markets in ecoinvent based on IAM output data.

    :ivar database: database dictionary from :attr:`.NewDatabase.database`
    :ivar iam_data: xarray that contains IAM data, from :attr:`.NewDatabase.rdc`
    :ivar model: can be 'remind' or 'image'. str from :attr:`.NewDatabase.model`
    :ivar pathway: file path to vehicle inventories
    :ivar year: year, from :attr:`.NewDatabase.year`
    :ivar version: ecoinvent database version
    :ivar relink: whether to relink supplier of datasets to better-fitted suppliers
    :ivar vehicle_type: "two-wheeler", "car", "bus" or "truck"
    :ivar has_fleet: whether `vehicle_type` has associated fleet data or not

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
        relink: bool,
        vehicle_type: str,
        has_fleet: bool,
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
            index,
        )
        self.version = version
        self.relink = relink
        self.vehicle_type = vehicle_type
        self.has_fleet = has_fleet
        self.database = database
        self.mapping = get_vehicles_mapping()

        if self.has_fleet:
            fleet_datasets = self.create_vehicle_markets()
            self.database.extend(fleet_datasets)
            self.add_to_index(fleet_datasets)

    def create_vehicle_markets(self) -> list:
        """
        Create vehicle market (fleet average) datasets.
        """

        # create and regionalize transport datasets
        filters = [
            ws.either(
                *[
                    ws.contains("name", size)
                    for size in self.mapping[self.vehicle_type]["sizes"]
                ]
            ),
            ws.equals("unit", self.mapping[self.vehicle_type]["unit"]),
        ]
        if self.vehicle_type != "truck":
            filters.append(ws.contains("name", self.vehicle_type))

        vehicle_datasets = list(ws.get_many(self.database, *filters))

        reduced_vehicle_datasets = list(
            set([(ds["name"], ds["reference product"]) for ds in vehicle_datasets])
        )

        for ds in reduced_vehicle_datasets:
            new_vehicle_datasets = self.fetch_proxies(
                subset=vehicle_datasets,
                name=ds[0],
                ref_prod=ds[1],
            )

            self.database.extend(new_vehicle_datasets.values())

        fleet_act = []

        arr = None
        if self.vehicle_type == "car":
            arr = self.iam_data.trsp_cars
        if self.vehicle_type == "truck":
            arr = self.iam_data.trsp_trucks
        if self.vehicle_type == "bus":
            arr = self.iam_data.trsp_buses

        if arr is None:
            return []

        # rename coordinates along the powertrian dimension
        rev_powertrain = {v: k for k, v in self.mapping["powertrain"].items()}
        arr.coords["powertrain"] = [
            rev_powertrain[p] for p in arr.coords["powertrain"].values
        ]

        vehicle_datasets = list(ws.get_many(self.database, *filters))

        fleet_act.extend(
            create_fleet_vehicles(
                datasets=vehicle_datasets,
                vehicle_type=self.vehicle_type,
                year=self.year,
                model=self.model,
                version=self.version,
                system_model=self.system_model,
                scenario=self.scenario,
                regions=self.regions,
                arr=arr,
            )
        )

        # if trucks, need to reconnect everything
        # loop through datasets that use truck transport
        if self.vehicle_type == "truck":
            list_created_trucks = [(a["name"], a["location"]) for a in fleet_act]
            for dataset in ws.get_many(
                self.database,
                ws.doesnt_contain_any("name", ["freight, lorry"]),
                ws.exclude(ws.equals("unit", "ton kilometer")),
            ):
                for exc in ws.technosphere(
                    dataset,
                    ws.contains("name", "transport, freight, lorry"),
                    ws.equals("unit", "ton kilometer"),
                ):
                    key = [
                        k
                        for k in self.mapping["truck"]["old_trucks"][self.model]
                        if k.lower() in exc["name"].lower()
                    ][0]

                    if "input" in exc:
                        del exc["input"]

                    if dataset["unit"] == "kilogram":
                        name = f"{self.mapping['truck']['old_trucks'][self.model][key]}, long haul"
                        loc = self.geo.ecoinvent_to_iam_location(dataset["location"])

                        if (name, loc) in list_created_trucks:
                            exc["name"] = name
                        else:
                            exc["name"] = (
                                f"transport, freight, lorry, unspecified, long haul"
                            )

                    else:
                        exc["name"] = (
                            "transport, freight, lorry, unspecified, long haul"
                        )

                    exc["product"] = "transport, freight, lorry"
                    exc["location"] = self.geo.ecoinvent_to_iam_location(
                        dataset["location"]
                    )

        return fleet_act
