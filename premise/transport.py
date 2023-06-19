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

from . import DATA_DIR, INVENTORY_DIR
from .inventory_imports import VariousVehicles
from .transformation import BaseTransformation, IAMDataCollection
from .utils import eidb_label

FILEPATH_FLEET_COMP = (
    DATA_DIR / "iam_output_files" / "fleet_files" / "fleet_all_vehicles.csv"
)
FILEPATH_IMAGE_TRUCKS_FLEET_COMP = (
    DATA_DIR / "iam_output_files" / "fleet_files" / "image_fleet_trucks.csv"
)
FILEPATH_TWO_WHEELERS = INVENTORY_DIR / "lci-two_wheelers.xlsx"
FILEPATH_TRUCKS = INVENTORY_DIR / "lci-trucks.xlsx"
FILEPATH_BUSES = INVENTORY_DIR / "lci-buses.xlsx"
FILEPATH_PASS_CARS = INVENTORY_DIR / "lci-pass_cars.xlsx"
FILEPATH_TRUCK_LOAD_FACTORS = DATA_DIR / "transport" / "avg_load_factors.yaml"
FILEPATH_VEHICLES_MAP = DATA_DIR / "transport" / "vehicles_map.yaml"


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
    print("Create fleet average vehicles...")

    vehicles_map = get_vehicles_mapping()

    if model == "remind":
        constr_year_map = {
            year: int(year.split("-")[-1])
            for year in arr.coords["construction_year"].values
        }
    else:
        constr_year_map = {
            year: year for year in arr.coords["construction_year"].values
        }

    # fleet data does not go below 2015
    if year < 2015:
        year = 2015
        print(
            "Vehicle fleet data is not available before 2015. "
            "Hence, 2015 is used as fleet year."
        )

    # fleet data does not go beyond 2050
    if year > 2050:
        year = 2050
        print(
            "Vehicle fleet data is not available beyond 2050. "
            "Hence, 2050 is used as fleet year."
        )

    # We filter electric vehicles by year of manufacture
    available_years = np.arange(2000, 2055, 5)
    ref_year = min(available_years, key=lambda x: abs(x - year))

    available_ds, d_names, cycle_type = [], {}, None

    for dataset in datasets:
        if dataset["name"].startswith("transport, "):
            if vehicle_type == "bus":
                if len(dataset["name"].split(", ")) == 6:
                    if "battery electric" in dataset["name"].split(", ")[2]:
                        _, _, pwt, _, size, year = dataset["name"].split(", ")

                    else:
                        _, _, pwt, size, year, _ = dataset["name"].split(", ")
                else:
                    _, _, pwt, size, year = dataset["name"].split(", ")

            elif vehicle_type == "truck":
                if len(dataset["name"].split(", ")) == 8:
                    if "battery electric" in dataset["name"].split(", ")[3]:
                        _, _, _, pwt, _, size, year, cycle_type = dataset["name"].split(
                            ", "
                        )

                    else:
                        _, _, _, pwt, size, year, _, cycle_type = dataset["name"].split(
                            ", "
                        )
                else:
                    _, _, _, pwt, size, year, cycle_type = dataset["name"].split(", ")

                size = size.replace(" gross weight", "")

            else:
                if len(dataset["name"].split(", ")) == 6:
                    if dataset["name"].split(", ")[2] == "battery electric":
                        _, _, pwt, _, size, year = dataset["name"].split(", ")
                    else:
                        _, _, pwt, size, year, _ = dataset["name"].split(", ")
                else:
                    _, _, pwt, size, year = dataset["name"].split(", ")

            if vehicle_type == "truck":
                d_names[
                    (vehicles_map["powertrain"][pwt], size, int(year), cycle_type)
                ] = (
                    dataset["name"],
                    dataset["reference product"],
                    dataset["unit"],
                )
                available_ds.append((vehicles_map["powertrain"][pwt], size, int(year)))

            else:
                d_names[(vehicles_map["powertrain"][pwt], size, int(year))] = (
                    dataset["name"],
                    dataset["reference product"],
                    dataset["unit"],
                )
                available_ds.append((vehicles_map["powertrain"][pwt], size, int(year)))

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
        if region not in arr.coords["region"].values:
            fleet_region = d_missing_regions[region]
        else:
            fleet_region = region

        sizes = [
            s
            for s in vehicles_map[vehicle_type]["sizes"]
            if s in arr.coords["size"].values
        ]

        sel = arr.sel(region=fleet_region, size=sizes, year=ref_year)

        total_km = sel.sum()

        if total_km > 0:
            if vehicle_type == "truck":
                driving_cycles = ["regional delivery", "long haul"]
            else:
                driving_cycles = [""]

            for driving_cycle in driving_cycles:
                name = (
                    f"{vehicles_map[vehicle_type]['name']}, unspecified, {driving_cycle}"
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
                    "database": eidb_label(
                        model, scenario, year, version, system_model
                    ),
                    "comment": f"Fleet-average vehicle for the year {year}, "
                    f"for the region {region}.",
                }

                for size in sizes:
                    for construction_year in sel.coords["construction_year"].values:
                        for pwt in sel.coords["powertrain"].values:
                            indiv_km = sel.sel(
                                size=size,
                                construction_year=construction_year,
                                powertrain=pwt,
                            )
                            if (
                                indiv_km > 0
                                and (pwt, size, constr_year_map[construction_year])
                                in available_ds
                            ):
                                indiv_share = (indiv_km / total_km).values.item(0)

                                if vehicle_type == "truck":
                                    load = avg_load[vehicle_type][driving_cycle][size]
                                    to_look_for = (
                                        pwt,
                                        size,
                                        constr_year_map[construction_year],
                                        driving_cycle,
                                    )

                                else:
                                    load = 1
                                    to_look_for = (
                                        pwt,
                                        size,
                                        constr_year_map[construction_year],
                                    )

                                if to_look_for in d_names:
                                    name, ref, unit = d_names[to_look_for]

                                    act["exchanges"].append(
                                        {
                                            "name": name,
                                            "product": ref,
                                            "unit": unit,
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
                        total_size_km = sel.sel(size=size).sum()

                        if total_size_km > 0:
                            name = (
                                f"{vehicles_map[vehicle_type]['name']}, {size} gross weight, "
                                f"unspecified powertrain, {driving_cycle}"
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
                                    model, scenario, year, version, system_model
                                ),
                                "comment": f"Fleet-average vehicle for the year {year}, for the region {region}.",
                            }

                            for construction_year in sel.coords[
                                "construction_year"
                            ].values:
                                for pwt in sel.coords["powertrain"].values:
                                    indiv_km = sel.sel(
                                        size=size,
                                        construction_year=construction_year,
                                        powertrain=pwt,
                                    )
                                    if (
                                        indiv_km > 0
                                        and (
                                            pwt,
                                            size,
                                            constr_year_map[construction_year],
                                        )
                                        in available_ds
                                    ):
                                        indiv_share = (
                                            indiv_km / total_size_km
                                        ).values.item(0)
                                        load = avg_load[vehicle_type][driving_cycle][
                                            size
                                        ]
                                        to_look_for = (
                                            pwt,
                                            size,
                                            constr_year_map[construction_year],
                                            driving_cycle,
                                        )
                                        if to_look_for in d_names:
                                            name, ref, unit = d_names[to_look_for]

                                            act["exchanges"].append(
                                                {
                                                    "name": name,
                                                    "product": ref,
                                                    "unit": unit,
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
        self.version = version
        self.relink = relink
        self.vehicle_type = vehicle_type
        self.has_fleet = has_fleet

    def generate_vehicles_datasets(self):
        """
        Generate vehicles datasets.
        """

        if self.vehicle_type == "car":
            filepath = FILEPATH_PASS_CARS
        elif self.vehicle_type == "truck":
            filepath = FILEPATH_TRUCKS
        elif self.vehicle_type == "bus":
            filepath = FILEPATH_BUSES
        else:
            filepath = FILEPATH_TWO_WHEELERS

        # load carculator inventories
        various_veh = VariousVehicles(
            database=self.database,
            version_in="3.7",
            version_out=self.version,
            path=filepath,
            year=self.year,
            regions=self.regions,
            model=self.model,
            scenario=self.scenario,
            vehicle_type=self.vehicle_type,
            has_fleet=True,
            system_model=self.system_model,
        )

        various_veh.prepare_inventory()

        return various_veh

    def create_vehicle_markets(self):
        """
        Create vehicle market (fleet average) datasets.
        """

        # create datasets
        datasets = self.generate_vehicles_datasets()

        list_vehicles = [
            "Bicycle,",
            "Kick-scooter,",
            "Moped,",
            "Scooter,",
            "Motorbike,",
            "urban delivery",
            "passenger bus",
        ]

        # We filter vehicles by year of manufacture
        available_years = [2020, 2030, 2040, 2050]
        closest_year = min(available_years, key=lambda x: abs(x - self.year))
        fleet_act = None

        if self.has_fleet:
            datasets.import_db.data = [
                dataset
                for dataset in datasets.import_db.data
                if not any(
                    z for z in re.findall(r"\d+", dataset["name"]) if int(z) > self.year
                )
                and not any(
                    z in dataset["name"]
                    for z in ["label-certified electricity", "label-certified gas"]
                )
            ]

            if self.vehicle_type == "car":
                arr = self.iam_data.trsp_cars
            elif self.vehicle_type == "truck":
                arr = self.iam_data.trsp_trucks
            elif self.vehicle_type == "bus":
                arr = self.iam_data.trsp_buses
            else:
                arr = None

            if arr is not None:
                fleet_act = create_fleet_vehicles(
                    datasets.import_db.data,
                    vehicle_type=self.vehicle_type,
                    year=self.year,
                    model=self.model,
                    version=self.version,
                    system_model=self.system_model,
                    scenario=self.scenario,
                    regions=self.regions,
                    arr=arr,
                )

                # cleaning up
                # we remove vehicles that
                # are not used by fleet vehicles
                fleet_vehicles = []

                for dataset in fleet_act:
                    for exchange in dataset["exchanges"]:
                        if exchange["type"] == "technosphere":
                            fleet_vehicles.append(exchange["name"])

                datasets.import_db.data = [
                    a
                    for a in datasets.import_db.data
                    if not a["name"].startswith("transport, ")
                    or a["name"] in fleet_vehicles
                ]

                datasets.import_db.data.extend(fleet_act)

        else:
            datasets.import_db.data = [
                dataset
                for dataset in datasets.import_db.data
                if not any(vehicle in dataset["name"] for vehicle in list_vehicles)
                or (
                    str(closest_year) in dataset["name"]
                    and "label-certified electricity" not in dataset["name"]
                )
            ]

            # remove the year in the name
            str_to_replace = ", " + str(closest_year)
            for dataset in datasets.import_db.data:
                if str_to_replace in dataset["name"]:
                    dataset["name"] = dataset["name"].replace(str_to_replace, "")
                    for exc in dataset["exchanges"]:
                        if str_to_replace in exc["name"]:
                            exc["name"] = exc["name"].replace(str_to_replace, "")

        list_new_ds = []

        # create regional variants
        for dataset in datasets.import_db.data:
            if (
                "transport, " in dataset["name"]
                and "unspecified" not in dataset["name"]
            ):
                for region in self.regions:
                    new_ds = wt.copy_to_new_location(dataset, region)

                    for exc in ws.production(new_ds):
                        if "input" in exc:
                            exc.pop("input")

                    if self.relink:
                        new_ds = self.relink_technosphere_exchanges(
                            new_ds,
                        )

                    list_new_ds.append(new_ds)

        datasets.import_db.data.extend(list_new_ds)

        # remove empty fields
        for dataset in datasets.import_db.data:
            for key, value in list(dataset.items()):
                if not value:
                    del dataset[key]

        # if trucks, need to reconnect everything
        # loop through datasets that use truck transport
        if self.vehicle_type == "truck":
            vehicles_map = get_vehicles_mapping()
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
                        for k in vehicles_map["truck"]["old_trucks"][self.model]
                        if k.lower() in exc["name"].lower()
                    ][0]

                    if "input" in exc:
                        del exc["input"]

                    if dataset["unit"] == "kilogram":
                        if exc["amount"] * 1000 <= 450:
                            name = f"{vehicles_map['truck']['old_trucks'][self.model][key]}, regional delivery"
                            cycle = ", regional delivery"
                        else:
                            name = f"{vehicles_map['truck']['old_trucks'][self.model][key]}, long haul"
                            cycle = ", long haul"

                        loc = self.geo.ecoinvent_to_iam_location(dataset["location"])
                        if (name, loc) in list_created_trucks:
                            exc["name"] = name

                        else:
                            exc["name"] = (
                                "transport, freight, lorry, unspecified" + cycle
                            )
                    else:
                        exc[
                            "name"
                        ] = "transport, freight, lorry, unspecified, long haul"

                    exc["product"] = "transport, freight, lorry"
                    exc["location"] = self.geo.ecoinvent_to_iam_location(
                        dataset["location"]
                    )

        self.database = datasets.merge_inventory()
