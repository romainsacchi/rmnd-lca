"""
transport.py contains the class Transport, which takes care of importing inventories
for a number of different vehicle types, and create fleet average vehicles based on
IAM data, and integrate them into the database.
"""

import uuid
from typing import Any, Dict, List, Union

import numpy as np
import xarray as xr
import yaml
from wurst import searching as ws

from .activity_maps import InventorySet
from .filesystem_constants import DATA_DIR, IAM_OUTPUT_DIR
from .logger import create_logger
from .transformation import BaseTransformation, IAMDataCollection
from .utils import eidb_label, rescale_exchanges
from .validation import CarValidation, TruckValidation

logger = create_logger("transport")

FILEPATH_FLEET_COMP = IAM_OUTPUT_DIR / "fleet_files" / "fleet_all_vehicles.csv"
FILEPATH_IMAGE_TRUCKS_FLEET_COMP = (
    IAM_OUTPUT_DIR / "fleet_files" / "image_fleet_trucks.csv"
)
FILEPATH_TRUCK_LOAD_FACTORS = DATA_DIR / "transport" / "avg_load_factors.yaml"
FILEPATH_VEHICLES_MAP = DATA_DIR / "transport" / "vehicles_map.yaml"


def _update_vehicles(scenario, vehicle_type, version, system_model):
    has_fleet = False
    if vehicle_type == "car":
        if hasattr(scenario["iam data"], "passenger_car_markets"):
            has_fleet = True
    elif vehicle_type == "truck":
        if hasattr(scenario["iam data"], "roadfreight_markets"):
            has_fleet = True
    elif vehicle_type == "train":
        if hasattr(scenario["iam data"], "railfreight_markets"):
            has_fleet = True
    elif vehicle_type == "bus":
        if hasattr(scenario["iam data"], "bus_markets"):
            has_fleet = True
    elif vehicle_type == "two-wheeler":
        if hasattr(scenario["iam data"], "two_wheelers_markets"):
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

    validation_func = {
        "car": CarValidation,
        "truck": TruckValidation,
    }

    if vehicle_type in validation_func:
        validate = validation_func[vehicle_type](
            model=scenario["model"],
            scenario=scenario["pathway"],
            year=scenario["year"],
            regions=scenario["iam data"].regions,
            database=trspt.database,
            iam_data=scenario["iam data"],
        )
        validate.run_checks()

    return scenario


def get_battery_size() -> dict:
    """
    Return a dictionary that contains the size of the battery
    for each vehicle type and powertrain.
    :return: dictionary with battery sizes
    """
    with open(
        DATA_DIR / "transport" / "battery_size.yaml", "r", encoding="utf-8"
    ) as stream:
        out = yaml.safe_load(stream)
        return out


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
    mapping: dict,
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

    for region in regions:
        if year in arr.coords["year"].values:
            region_size_fleet = arr.sel(region=region, year=year)

        else:
            region_size_fleet = arr.sel(region=region).interp(year=year)

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

            for vehicle in arr.coords["variables"].values:
                indiv_km = region_size_fleet.sel(
                    variables=vehicle,
                )
                if indiv_km > 0:
                    indiv_share = (indiv_km / total_km).values.item(0)

                    try:
                        name = mapping[vehicle]
                    except KeyError:
                        print(mapping)
                        name = vehicle

                    if isinstance(name, set):
                        # check if length of set is 1
                        if len(name) == 1:
                            # if so, take the only element
                            name = next(iter(name))

                    try:
                        vehicle_dataset = ws.get_one(
                            datasets,
                            ws.equals("name", name),
                            ws.equals("location", region),
                        )
                    except ws.NoResults:
                        print(f"Could not find dataset for {name} in {region}.")
                        continue
                    except ws.MultipleResults:
                        print(f"Multiple datasets found for {name} in {region}.")
                        continue

                    act["exchanges"].append(
                        {
                            "name": vehicle_dataset["name"],
                            "product": vehicle_dataset["reference product"],
                            "unit": vehicle_dataset["unit"],
                            "location": region,
                            "type": "technosphere",
                            "amount": indiv_share,
                        }
                    )

            if len(act["exchanges"]) > 1:
                list_act.append(act)

            # also create size-specific fleet vehicles
            if vehicle_type == "truck":
                sizes = ["3.5t", "7.5t", "18t", "26t", "40t"]
                for size in sizes:
                    total_size_km = region_size_fleet.sel(
                        variables=[
                            v for v in arr.coords["variables"].values if size in v
                        ]
                    ).sum()

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

                        for pwt in [
                            v for v in arr.coords["variables"].values if size in v
                        ]:
                            indiv_km = region_size_fleet.sel(
                                variables=pwt,
                            )
                            if indiv_km > 0:
                                indiv_share = (indiv_km / total_size_km).values.item(0)

                                name = mapping[pwt]
                                if isinstance(name, set):
                                    # check if length of set is 1
                                    if len(name) == 1:
                                        # if so, take the only element
                                        name = next(iter(name))

                                vehicle_dataset = ws.get_one(
                                    datasets,
                                    ws.equals("name", name),
                                    ws.equals("location", region),
                                )

                                act["exchanges"].append(
                                    {
                                        "name": vehicle_dataset["name"],
                                        "product": vehicle_dataset["reference product"],
                                        "unit": vehicle_dataset["unit"],
                                        "location": region,
                                        "type": "technosphere",
                                        "amount": indiv_share,
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

        mapping = InventorySet(database=database, version=version, model=model)
        self.vehicle_map = mapping.generate_transport_map(transport_type=vehicle_type)
        self.rev_map = {next(iter(v)): k for k, v in self.vehicle_map.items()}
        self.vehicle_fuel_map = mapping.generate_vehicle_fuel_map(
            transport_type=vehicle_type
        )
        self.battery_size = get_battery_size()

        # check if vehicle map is empty
        for v in self.vehicle_map.values():
            if not v:
                print(f"Vehicle map is empty for {self.vehicle_type}.")

        if self.has_fleet:
            fleet_datasets = self.create_vehicle_markets()
            self.database.extend(fleet_datasets)
            self.add_to_index(fleet_datasets)

    def create_vehicle_markets(self) -> list:
        """
        Create vehicle market (fleet average) datasets.
        """

        # create and regionalize transport datasets
        vehicle_datasets = list(
            ws.get_many(
                self.database,
                ws.either(
                    *[
                        ws.equals("name", v)
                        for name in self.vehicle_map.values()
                        for v in name
                    ]
                ),
            )
        )

        new_datasets = []

        for ds in list(
            set([(v["name"], v["reference product"]) for v in vehicle_datasets])
        ):
            new_datasets.extend(
                self.fetch_proxies(
                    subset=vehicle_datasets,
                    name=ds[0],
                    ref_prod=ds[1],
                ).values()
            )

        for new_ds in new_datasets:
            new_ds = self.adjust_transport_efficiency(new_ds)

            if not self.is_in_index(new_ds):
                self.add_to_index(new_ds)
                self.database.append(new_ds)

            else:
                print(
                    f"Dataset {new_ds['name'], new_ds['location']} already in the database."
                )

        fleet_act = []

        arr = None
        if self.vehicle_type == "two-wheeler":
            arr = self.iam_data.two_wheelers_markets
        if self.vehicle_type == "car":
            arr = self.iam_data.passenger_car_markets
        if self.vehicle_type == "truck":
            arr = self.iam_data.roadfreight_markets
        if self.vehicle_type == "bus":
            arr = self.iam_data.bus_markets
        if self.vehicle_type == "train":
            arr = self.iam_data.railfreight_markets

        if arr is None:
            return []

        fleet_act.extend(
            create_fleet_vehicles(
                datasets=new_datasets,
                vehicle_type=self.vehicle_type,
                year=self.year,
                model=self.model,
                version=self.version,
                system_model=self.system_model,
                scenario=self.scenario,
                regions=self.regions,
                arr=arr,
                mapping=self.vehicle_map,
            )
        )

        # if trucks, adjust battery size
        if self.vehicle_type == "truck":
            for ds in ws.get_many(
                self.database,
                ws.contains("name", "battery electric"),
                ws.contains("name", "truck"),
                ws.equals("unit", "unit"),
            ):
                self.adjust_battery_size(ds)

        # if trucks, need to reconnect everything
        # loop through datasets that use truck transport
        if self.vehicle_type == "truck":

            list_created_vehicles = [(v["name"], v["location"]) for v in fleet_act]

            for dataset in ws.get_many(
                self.database,
                ws.doesnt_contain_any("name", ["freight, lorry"]),
                ws.exclude(ws.equals("unit", "ton kilometer")),
            ):
                for exc in ws.technosphere(
                    dataset,
                    ws.either(
                        *[
                            ws.equals("name", v)
                            for v in self.mapping["truck"]["old_trucks"]
                        ]
                    ),
                    ws.equals("unit", "ton kilometer"),
                ):

                    new_name = self.mapping["truck"]["old_trucks"][exc["name"]][
                        self.model
                    ]
                    new_loc = self.geo.ecoinvent_to_iam_location(dataset["location"])

                    if (new_name, new_loc) in list_created_vehicles:
                        exc["name"] = new_name
                        exc["product"] = "transport, freight, lorry"
                        exc["location"] = new_loc
                    else:
                        print(f"Could not find dataset for {new_name} in {new_loc}.")
                        exc["name"] = (
                            "transport, freight, lorry, unspecified, long haul"
                        )
                        exc["product"] = "transport, freight, lorry"
                        exc["location"] = "World"

            # also we need to empty the old transport datasets
            for dataset in ws.get_many(
                self.database,
                ws.either(
                    *[ws.equals("name", v) for v in self.mapping["truck"]["old_trucks"]]
                ),
            ):
                dataset["exchanges"] = [
                    e for e in dataset["exchanges"] if e["type"] == "production"
                ]
                dataset["comment"] = (
                    "This dataset has been replaced by new fleet-average vehicles."
                )

                # add new truck as exchange
                new_name = self.mapping["truck"]["old_trucks"][dataset["name"]][
                    self.model
                ]
                new_loc = self.geo.ecoinvent_to_iam_location(dataset["location"])

                if (new_name, new_loc) in list_created_vehicles:
                    new_exc = {
                        "name": new_name,
                        "product": "transport, freight, lorry",
                        "unit": "ton kilometer",
                        "location": new_loc,
                        "type": "technosphere",
                        "amount": 1,
                        "uncertainty type": 0,
                    }
                else:
                    print(f"Could not find dataset for {new_name} in {new_loc}.")
                    new_exc = {
                        "name": "transport, freight, lorry, unspecified, long haul",
                        "product": "transport, freight, lorry",
                        "unit": "ton kilometer",
                        "location": "World",
                        "type": "technosphere",
                        "amount": 1,
                        "uncertainty type": 0,
                    }
                dataset["exchanges"].append(new_exc)

        return fleet_act

    def adjust_transport_efficiency(self, dataset):
        """
        Adjust transport efficiency of transport datasets based on IAM data.

        :param dataset: dataset to adjust
        :return: dataset with adjusted transport efficiency
        """

        if self.vehicle_type == "car":
            data = self.iam_data.passenger_car_efficiencies
        elif self.vehicle_type == "truck":
            data = self.iam_data.roadfreight_efficiencies
        elif self.vehicle_type == "bus":
            data = self.iam_data.bus_efficiencies
        elif self.vehicle_type == "train":
            data = self.iam_data.railfreight_efficiencies
        elif self.vehicle_type == "two-wheeler":
            data = self.iam_data.two_wheelers_efficiencies
        else:
            raise ValueError("Unknown vehicle type.")

        if data is None:
            return dataset

        variable = self.rev_map[dataset["name"]]

        scaling_factor = 1 / self.find_iam_efficiency_change(
            data=data,
            variable=variable,
            location=dataset["location"],
        )

        if scaling_factor != 1:
            dataset = rescale_exchanges(
                dataset,
                scaling_factor,
                technosphere_filters=[
                    ws.either(
                        *[
                            ws.contains("name", v)
                            for v in self.vehicle_fuel_map[variable]
                        ]
                    )
                ],
            )
            if "log parameters" not in dataset:
                dataset["log parameters"] = {}

            dataset["log parameters"].update({"efficiency change": scaling_factor})

        self.write_log(dataset)

        return dataset

    def adjust_battery_size(self, ds):
        """
        Adjust battery size for truck datasets.
        """

        # detect size in name
        size = [s for s in self.battery_size["truck"] if s in ds["name"]][0]

        if self.year <= min(self.battery_size["truck"][size].keys()):
            mean_battery_size = self.battery_size["truck"][size][
                min(self.battery_size["truck"][size].keys())
            ]["mean"]
            min_battery_size = self.battery_size["truck"][size][
                min(self.battery_size["truck"][size].keys())
            ]["min"]
            max_battery_size = self.battery_size["truck"][size][
                min(self.battery_size["truck"][size].keys())
            ]["max"]
        elif self.year >= max(self.battery_size["truck"][size].keys()):
            mean_battery_size = self.battery_size["truck"][size][
                max(self.battery_size["truck"][size].keys())
            ]["mean"]
            min_battery_size = self.battery_size["truck"][size][
                max(self.battery_size["truck"][size].keys())
            ]["min"]
            max_battery_size = self.battery_size["truck"][size][
                max(self.battery_size["truck"][size].keys())
            ]["max"]
        else:
            mean_battery_size = np.interp(
                self.year,
                list(self.battery_size["truck"][size].keys()),
                [v["mean"] for v in self.battery_size["truck"][size].values()],
            )
            min_battery_size = np.interp(
                self.year,
                list(self.battery_size["truck"][size].keys()),
                [v["min"] for v in self.battery_size["truck"][size].values()],
            )
            max_battery_size = np.interp(
                self.year,
                list(self.battery_size["truck"][size].keys()),
                [v["max"] for v in self.battery_size["truck"][size].values()],
            )

        for exc in ws.technosphere(ds, ws.contains("name", "market for battery")):
            exc["amount"] = mean_battery_size
            exc["uncertainty type"] = 5
            exc["loc"] = exc["amount"]
            exc["minimum"] = min_battery_size
            exc["maximum"] = max_battery_size

        ds["comment"] += f" Battery size adjusted to {mean_battery_size} kWh."

    def write_log(self, dataset, status="created"):
        """
        Write log file.
        """

        logger.info(
            f"{status}|{self.model}|{self.scenario}|{self.year}|"
            f"{dataset['name']}|{dataset['location']}|"
            f"{dataset.get('log parameters', {}).get('efficiency change', '')}"
        )
