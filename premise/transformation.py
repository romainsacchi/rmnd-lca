"""
transformation.py contains the base class TransformationBase,
used by other classes (e.g. Transport, Electricity, Steel, Cement, etc.).
It provides basic methods usually used for electricity, cement, steel sector transformation
on the wurst database.
"""

import logging.config
import uuid
from collections import defaultdict
from copy import deepcopy
from functools import lru_cache
from itertools import groupby, product
from pathlib import Path
from typing import Any, Dict, List, Set, Tuple, Union

import numpy as np
import xarray as xr
import yaml
from _operator import itemgetter
from constructive_geometries import resolved_row
from wurst import reference_product, rescale_exchange
from wurst import searching as ws
from wurst import transformations as wt

from .activity_maps import InventorySet
from .data_collection import IAMDataCollection
from .geomap import Geomap
from .utils import DATA_DIR, get_fuel_properties

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

logger = logging.getLogger("module")


def get_suppliers_of_a_region(
    database: List[dict],
    locations: List[str],
    names: List[str],
    reference_prod: str,
    unit: str,
    exclude: List[str] = None,
) -> filter:
    """
    Return a list of datasets, for which the location, name,
    reference product and unit correspond to the region and name
    given, respectively.

    :param database: database to search
    :param locations: list of locations
    :param names: names of datasets
    :param unit: unit of dataset
    :param reference_prod: reference product of dataset
    :return: list of wurst datasets
    :param exclude: list of terms to exclude
    """

    filters = [
        ws.either(*[ws.contains("name", supplier) for supplier in names]),
        ws.either(*[ws.equals("location", loc) for loc in locations]),
        ws.contains("reference product", reference_prod),
        ws.equals("unit", unit),
    ]

    if exclude:
        filters.append(ws.doesnt_contain_any("name", exclude))

    return ws.get_many(
        database,
        *filters,
    )


def get_shares_from_production_volume(
    ds_list: Union[Dict[str, Any], List[Dict[str, Any]]],
) -> Dict[Tuple[Any, Any, Any, Any], float]:
    """
    Return shares of supply of each dataset in `ds_list`
    based on respective production volumes
    :param ds_list: list of datasets
    :return: dictionary with (dataset name, dataset location, ref prod, unit) as keys, shares as values. Shares total 1.
    """

    if not isinstance(ds_list, list):
        ds_list = [ds_list]

    dict_act = {}
    total_production_volume = 0

    for act in ds_list:
        production_volume = 0

        if "production volume" in act:
            production_volume = max(float(act["production volume"]), 1e-9)
        else:
            for exc in ws.production(act):
                # even if non-existent, we set a minimum value of 1e-9
                # because if not, we risk dividing by zero!!!
                production_volume = max(float(exc.get("production volume", 1e-9)), 1e-9)

        dict_act[
            (
                act["name"],
                act["location"],
                act["reference product"],
                act["unit"],
            )
        ] = production_volume
        total_production_volume += production_volume

    for dataset in dict_act:
        dict_act[dataset] /= total_production_volume

    return dict_act


def get_tuples_from_database(database: List[dict]) -> List[Tuple[str, str, str]]:
    """
    Return a list of tuples (name, reference product, location)
    for each dataset in database.
    :param database: wurst database
    :return: a list of tuples
    """
    return [
        (dataset["name"], dataset["reference product"], dataset["location"])
        for dataset in database
        if "has_downstream_consumer" not in dataset
    ]


def remove_exchanges(datasets_dict: Dict[str, dict], list_exc: List) -> Dict[str, dict]:
    """
    Returns the same `datasets_dict`, where the list of exchanges in these datasets
    has been filtered out: unwanted exchanges has been removed.

    :param datasets_dict: a dictionary with IAM regions as keys, datasets as value
    :param list_exc: list of names (e.g., ["coal", "lignite"]) which are checked against exchanges' names in the dataset
    :return: returns `datasets_dict` without the exchanges whose names check with `list_exc`
    """
    keep = lambda x: {
        key: value
        for key, value in x.items()
        if not any(ele in x.get("product", []) for ele in list_exc)
    }

    for region in datasets_dict:
        datasets_dict[region]["exchanges"] = [
            keep(exc) for exc in datasets_dict[region]["exchanges"]
        ]

    return datasets_dict


def new_exchange(exc, location, factor):
    copied_exc = deepcopy(exc)
    copied_exc["location"] = location
    return rescale_exchange(copied_exc, factor)


def get_gis_match(
    dataset,
    location,
    possible_locations,
    geomatcher,
    contained,
    exclusive,
    biggest_first,
):
    with resolved_row(possible_locations, geomatcher.geo) as g:
        func = g.contained if contained else g.intersects

        if dataset["location"] not in geomatcher.iam_regions:
            gis_match = func(
                location,
                include_self=True,
                exclusive=exclusive,
                biggest_first=biggest_first,
                only=possible_locations,
            )

        else:
            gis_match = geomatcher.iam_to_ecoinvent_location(dataset["location"])

    return gis_match


def allocate_inputs(exc, lst):
    """
    Allocate the input exchanges in ``lst`` to ``exc``,
    using production volumes where possible, and equal splitting otherwise.
    Always uses equal splitting if ``RoW`` is present.
    """
    pvs = [reference_product(o).get("production volume") or 0 for o in lst]

    if all((x > 0 for x in pvs)):
        # Allocate using production volume
        total = sum(pvs)
    else:
        # Allocate evenly
        total = len(lst)
        pvs = [1 for _ in range(total)]

    return [
        new_exchange(exc, obj["location"], factor / total)
        for obj, factor in zip(lst, pvs)
    ], [p / total for p in pvs]


def filter_out_results(
    item_to_look_for: str, results: List[dict], field_to_look_at: str
) -> List[dict]:
    """Filters a list of results by a given field"""
    return [r for r in results if item_to_look_for not in r[field_to_look_at]]


class BaseTransformation:
    """
    Base transformation class.

    :ivar database: wurst database
    :ivar iam_data: IAMDataCollection object_
    :ivar model: IAM model
    :ivar year: database year
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
        cache: dict = None,
    ) -> None:
        self.database: List[dict] = database
        self.iam_data: IAMDataCollection = iam_data
        self.model: str = model
        self.regions: List[str] = iam_data.regions
        self.geo: Geomap = Geomap(model=model)
        self.scenario: str = pathway
        self.year: int = year
        self.version: str = version
        self.fuels_specs: dict = get_fuel_properties()
        mapping = InventorySet(self.database)
        self.cement_fuels_map: dict = mapping.generate_cement_fuels_map()
        self.fuel_map: Dict[str, Set] = mapping.generate_fuel_map()
        self.system_model: str = system_model
        self.modified_datasets = modified_datasets
        self.cache: dict = cache or {}

        # reverse the fuel map to get a mapping from ecoinvent to premise
        self.fuel_map_reverse: Dict = {}

        for key, value in self.fuel_map.items():
            for v in list(value):
                self.fuel_map_reverse[v] = key

        self.material_map: Dict[str, Set] = mapping.generate_material_map()
        self.ecoinvent_to_iam_loc: Dict[str, str] = {
            loc: self.geo.ecoinvent_to_iam_location(loc)
            for loc in self.get_ecoinvent_locs()
        }

    @lru_cache
    def select_multiple_suppliers(
        self,
        possible_names,
        dataset_location,
        look_for=None,
        blacklist=None,
    ):
        """
        Select multiple suppliers for a specific fuel.
        """

        # We have several potential fuel suppliers
        # We will look up their respective production volumes
        # And include them proportionally to it

        ecoinvent_regions = self.geo.iam_to_ecoinvent_location(dataset_location)

        possible_locations = [
            dataset_location,
            [*ecoinvent_regions],
            "RoW",
            "Europe without Switzerland",
            "RER",
            "GLO",
        ]

        suppliers, counter = [], 0

        extra_filters = []
        if look_for:
            extra_filters.append(
                ws.either(
                    *[ws.contains("reference product", item) for item in look_for]
                )
            )
        if blacklist:
            extra_filters.append(
                ws.exclude(
                    ws.either(*[ws.contains("reference product", x) for x in blacklist])
                )
            )
            extra_filters.append(
                ws.exclude(ws.either(*[ws.contains("name", x) for x in blacklist]))
            )

        try:
            while not suppliers:
                suppliers = list(
                    ws.get_many(
                        self.database,
                        ws.either(
                            *[ws.contains("name", sup) for sup in possible_names]
                        ),
                        ws.either(
                            *[
                                ws.equals("location", item)
                                for item in possible_locations[counter]
                            ]
                        )
                        if isinstance(possible_locations[counter], list)
                        else ws.equals("location", possible_locations[counter]),
                        *extra_filters,
                    )
                )
                counter += 1
        except IndexError:
            raise IndexError(
                "No supplier found for {} in {}, "
                "looking for terms: {} "
                "and with blacklist: {}".format(
                    possible_names, possible_locations, look_for, blacklist
                )
            )

        suppliers = get_shares_from_production_volume(suppliers)

        return suppliers

    def get_ecoinvent_locs(self) -> List[str]:
        """
        Rerun a list of unique locations in ecoinvent

        :return: list of locations
        :rtype: list
        """

        return list(set([a["location"] for a in self.database]))

    def update_ecoinvent_efficiency_parameter(
        self, dataset: dict, old_ei_eff: float, new_eff: float
    ) -> None:
        """
        Update the old efficiency value in the ecoinvent dataset by the newly calculated one.
        :param dataset: dataset
        :param old_ei_eff: conversion efficiency of the original ecoinvent dataset
        :param new_eff: new conversion efficiency
        :return: nothing. Modifies the `comment` and `parameters` fields of the dataset.
        """
        parameters = dataset["parameters"]
        possibles = ["efficiency", "efficiency_oil_country", "efficiency_electrical"]

        if any(i in dataset for i in possibles):
            for key in possibles:
                if key in parameters:
                    dataset["parameters"][key] = new_eff
        else:
            dataset["parameters"]["efficiency"] = new_eff

        if dataset["location"] in self.regions:
            iam_region = dataset["location"]
        else:
            iam_region = self.ecoinvent_to_iam_loc[dataset["location"]]

        new_txt = (
            f" 'premise' has modified the efficiency of this dataset, from an original "
            f"{int(old_ei_eff * 100)}% to {int(new_eff * 100)}%, according to IAM model {self.model.upper()}, scenario {self.scenario} "
            f"for the region {iam_region}."
        )

        if "comment" in dataset:
            dataset["comment"] += new_txt
        else:
            dataset["comment"] = new_txt

    def calculate_input_energy(
        self, fuel_name: str, fuel_amount: float, fuel_unit: str
    ) -> float:
        """
        Returns the amount of energy entering the conversion process, in MJ
        :param fuel_name: name of the liquid, gaseous or solid fuel
        :param fuel_amount: amount of fuel input
        :param fuel_unit: unit of fuel
        :return: amount of fuel energy, in MJ
        """

        # if fuel input other than MJ
        if fuel_unit in ["kilogram", "cubic meter", "kilowatt hour"]:
            try:
                lhv = self.fuels_specs[self.fuel_map_reverse[fuel_name]]["lhv"]

            except KeyError:
                lhv = 0
        else:
            lhv = 1

        # if already in MJ
        return fuel_amount * lhv

    def find_fuel_efficiency(
        self, dataset: dict, fuel_filters: List[str], energy_out: float
    ) -> float:
        """
        This method calculates the efficiency value set initially, in case it is not specified in the parameter
        field of the dataset. In Carma datasets, fuel inputs are expressed in megajoules instead of kilograms.

        :param dataset: a wurst dataset of an electricity-producing technology
        :param fuel_filters: wurst filter to filter fuel input exchanges
        :param energy_out: the amount of energy expect as output, in MJ
        :return: the efficiency value set initially
        """

        not_allowed = ["thermal"]
        key = []
        if "parameters" in dataset:
            key = list(
                key
                for key in dataset["parameters"]
                if "efficiency" in key and not any(item in key for item in not_allowed)
            )

        if len(key) > 0:
            if "parameters" in dataset:
                dataset["parameters"]["efficiency"] = dataset["parameters"][key[0]]
            else:
                dataset["parameters"] = {"efficiency": dataset["parameters"][key[0]]}

            return dataset["parameters"][key[0]]

        energy_input = np.sum(
            np.sum(
                np.asarray(
                    [
                        self.calculate_input_energy(
                            exc["name"], exc["amount"], exc["unit"]
                        )
                        for exc in dataset["exchanges"]
                        if exc["name"] in fuel_filters and exc["type"] == "technosphere"
                    ]
                )
            )
        )

        if energy_input != 0 and float(energy_out) != 0:
            current_efficiency = float(energy_out) / energy_input
        else:
            current_efficiency = np.nan

        if current_efficiency in (np.nan, np.inf):
            current_efficiency = 1

        if "parameters" in dataset:
            dataset["parameters"]["efficiency"] = current_efficiency
        else:
            dataset["parameters"] = {"efficiency": current_efficiency}

        return current_efficiency

    def get_iam_mapping(
        self, activity_map: dict, fuels_map: dict, technologies: list
    ) -> Dict[str, Any]:
        """
        Define filter functions that decide which wurst datasets to modify.
        :param activity_map: a dictionary that contains 'technologies' as keys and activity names as values.
        :param fuels_map: a dictionary that contains 'technologies' as keys and fuel names as values.
        :param technologies: a list of IAM technologies.
        :return: dictionary that contains filters and functions
        :rtype: dict
        """

        return {
            tech: {
                "IAM_eff_func": self.find_iam_efficiency_change,
                "current_eff_func": self.find_fuel_efficiency,
                "technology filters": activity_map[tech],
                "fuel filters": fuels_map[tech],
            }
            for tech in technologies
        }

    def region_to_proxy_dataset_mapping(
        self, name: str, ref_prod: str, regions: List[str] = None
    ) -> Dict[str, str]:
        d_map = {
            self.ecoinvent_to_iam_loc[d["location"]]: d["location"]
            for d in ws.get_many(
                self.database,
                ws.equals("name", name),
                ws.contains("reference product", ref_prod),
            )
            if d["location"] not in self.regions
        }

        if not regions:
            regions = self.regions

        if "RoW" in d_map.values():
            fallback_loc = "RoW"
        else:
            if "GLO" in d_map.values():
                fallback_loc = "GLO"
            else:
                fallback_loc = list(d_map.values())[0]

        return {region: d_map.get(region, fallback_loc) for region in regions}

    def fetch_proxies(
        self,
        name,
        ref_prod,
        production_variable=None,
        relink=True,
        regions=None,
        delete_original_dataset=False,
    ) -> Dict[str, dict]:
        """
        Fetch dataset proxies, given a dataset `name` and `reference product`.
        Store a copy for each IAM region.
        If a fitting ecoinvent location cannot be found for a given IAM region,
        fetch a dataset with a "RoW" location.
        Delete original datasets from the database.

        :param name: name of the datasets to find
        :param ref_prod: reference product of the datasets to find
        :param production_variable: name of variable in IAM data that refers to production volume
        :param relink: if `relink`, exchanges from the datasets will be relinked to
        the most geographically-appropriate providers from the database. This is computer-intensive.
        :param regions: regions to create proxy datasets for. if None, all regions are considered.
        :return: dictionary with IAM regions as keys, proxy datasets as values.
        """

        d_iam_to_eco = self.region_to_proxy_dataset_mapping(
            name=name, ref_prod=ref_prod, regions=regions
        )

        d_act = {}

        ds_name, ds_ref_prod = [None, None]

        for region in d_iam_to_eco:
            try:
                dataset = ws.get_one(
                    self.database,
                    ws.equals("name", name),
                    ws.contains("reference product", ref_prod),
                    ws.equals("location", d_iam_to_eco[region]),
                )
            except ws.MultipleResults as err:
                print(
                    err,
                    "A single dataset was expected, "
                    f"but found more than one for: {name, ref_prod}",
                )

            if (name, ref_prod, region, dataset["unit"]) not in self.modified_datasets[
                (self.model, self.scenario, self.year)
            ]["created"]:
                d_act[region] = wt.copy_to_new_location(dataset, region)
                d_act[region]["code"] = str(uuid.uuid4().hex)

                for exc in ws.production(d_act[region]):
                    if "input" in exc:
                        del exc["input"]

                if "input" in d_act[region]:
                    del d_act[region]["input"]

                if production_variable:
                    # Add `production volume` field
                    if isinstance(production_variable, str):
                        production_variable = [production_variable]

                    if all(
                        i in self.iam_data.production_volumes.variables
                        for i in production_variable
                    ):
                        prod_vol = (
                            self.iam_data.production_volumes.sel(
                                region=region, variables=production_variable
                            )
                            .interp(year=self.year)
                            .sum(dim="variables")
                            .values.item(0)
                        )

                    else:
                        prod_vol = 1
                else:
                    prod_vol = 1

                for prod in ws.production(d_act[region]):
                    prod["location"] = region
                    prod["production volume"] = prod_vol

                if relink:
                    d_act[region] = self.relink_technosphere_exchanges(d_act[region])

                ds_name = d_act[region]["name"]
                ds_ref_prod = d_act[region]["reference product"]

        # add dataset to emptied datasets list
        for ds in self.database:
            if (ds["name"] == ds_name) and (ds["reference product"] == ds_ref_prod):
                self.modified_datasets[(self.model, self.scenario, self.year)][
                    "emptied"
                ].append(
                    (
                        ds["name"],
                        ds["reference product"],
                        ds["location"],
                        ds["unit"],
                    )
                )

        # empty original datasets
        # and make them link to new regional datasets
        self.empty_original_datasets(
            name=ds_name,
            ref_prod=ds_ref_prod,
            loc_map=d_iam_to_eco,
            production_variable=production_variable,
            regions=regions,
        )

        if delete_original_dataset:
            # remove the dataset from `self.database`
            self.database = [
                ds
                for ds in self.database
                if not (
                    ds["name"] == ds_name and ds["reference product"] == ds_ref_prod
                )
            ]

        return d_act

    def empty_original_datasets(
        self,
        name: str,
        ref_prod: str,
        production_variable: str,
        loc_map: Dict[str, str],
        regions: List[str] = None,
    ) -> None:
        """
        Empty original ecoinvent dataset and introduce an input to the regional IAM
        dataset that geographically comprises it.
        :param name: dataset name
        :param ref_prod: dataset reference product
        :param loc_map: ecoinvent location to IAM location mapping for this activity
        :param production_variable: IAM production variable
        :return: Does not return anything. Just empties the original dataset.
        """

        regions = regions or self.regions
        if loc_map:
            mapping = defaultdict(set)
            for k, v in loc_map.items():
                if self.geo.ecoinvent_to_iam_location(v) in loc_map.keys():
                    mapping[v].add(self.geo.ecoinvent_to_iam_location(v))

        existing_datasets = ws.get_many(
            self.database,
            ws.equals("name", name),
            ws.contains("reference product", ref_prod),
            ws.doesnt_contain_any("location", regions),
        )

        for existing_ds in existing_datasets:
            if existing_ds["location"] in mapping:
                iam_locs = list(mapping[existing_ds["location"]])
            else:
                iam_locs = [self.ecoinvent_to_iam_loc[existing_ds["location"]]]
                iam_locs = [loc for loc in iam_locs if loc in regions]

            if iam_locs == ["World"]:
                iam_locs = [r for r in regions if r != "World"]

            # add tag
            existing_ds["has_downstream_consumer"] = False
            existing_ds["exchanges"] = [
                e for e in existing_ds["exchanges"] if e["type"] == "production"
            ]

            # for cases where external scenarios are used
            if "adjust efficiency" in existing_ds:
                del existing_ds["adjust efficiency"]

            if len(existing_ds["exchanges"]) == 0:
                print(
                    f"ISSUE: no exchanges found in {existing_ds['name']} "
                    f"in {existing_ds['location']}"
                )

            _ = lambda x: x if x != 0.0 else 1.0

            if len(iam_locs) == 1:
                existing_ds["exchanges"].append(
                    {
                        "name": existing_ds["name"],
                        "product": existing_ds["reference product"],
                        "amount": 1.0,
                        "unit": existing_ds["unit"],
                        "uncertainty type": 0,
                        "location": iam_locs[0],
                        "type": "technosphere",
                    }
                )
            else:
                for iam_loc in iam_locs:
                    if production_variable and all(
                        i in self.iam_data.production_volumes.variables.values.tolist()
                        for i in production_variable
                    ):
                        share = (
                            self.iam_data.production_volumes.sel(
                                region=iam_loc, variables=production_variable
                            )
                            .interp(year=self.year)
                            .sum(dim="variables")
                            .values.item(0)
                        ) / _(
                            self.iam_data.production_volumes.sel(
                                region=iam_locs, variables=production_variable
                            )
                            .interp(year=self.year)
                            .sum(dim=["variables", "region"])
                            .values.item(0)
                        )

                    else:
                        share = 1 / len(iam_locs)

                    if share > 0:
                        existing_ds["exchanges"].append(
                            {
                                "name": existing_ds["name"],
                                "product": existing_ds["reference product"],
                                "amount": share,
                                "unit": existing_ds["unit"],
                                "uncertainty type": 0,
                                "location": iam_loc,
                                "type": "technosphere",
                            }
                        )

            # add dataset to emptied datasets list
            self.modified_datasets[(self.model, self.scenario, self.year)][
                "emptied"
            ].append(
                (
                    existing_ds["name"],
                    existing_ds["reference product"],
                    existing_ds["location"],
                    existing_ds["unit"],
                )
            )

            # add log
            self.write_log(dataset=existing_ds, status="empty")

    def relink_datasets(
        self, excludes_datasets: List[str] = None, alt_names: List[str] = None
    ) -> None:
        """
        For a given exchange name, product and unit, change its location to an IAM location,
        to effectively link to the newly built market(s)/activity(ies).
        :param excludes_datasets: list of datasets to exclude from relinking
        :param alt_names: list of alternative names to use for relinking
        """

        if alt_names is None:
            alt_names = []
        if excludes_datasets is None:
            excludes_datasets = []

        # loop through the database
        # ignore datasets which name contains `name`
        for act in ws.get_many(
            self.database,
            ws.doesnt_contain_any("name", excludes_datasets),
        ):
            # and find exchanges of datasets to relink
            excs_to_relink = [
                exchange
                for exchange in act["exchanges"]
                if exchange["type"] == "technosphere"
                and (
                    exchange["name"],
                    exchange["product"],
                    exchange["location"],
                    exchange["unit"],
                )
                in self.modified_datasets[(self.model, self.scenario, self.year)][
                    "emptied"
                ]
            ]

            unique_excs_to_relink = set(
                (exc["name"], exc["product"], exc["location"], exc["unit"])
                for exc in excs_to_relink
            )

            new_exchanges = []

            for exc in unique_excs_to_relink:
                entry = None

                # check if already in cache
                if exc in self.cache.get(act["location"], {}).get(self.model, {}):
                    entry = self.cache[act["location"]][self.model][exc]

                    if isinstance(entry, tuple):
                        entry = [entry + (1.0,)]

                # not in cache, so find new candidates
                else:
                    names_to_look_for = [exc[0], *alt_names]
                    if exc[0].startswith("market group for"):
                        names_to_look_for.append(
                            exc[0].replace("market group for", "market for")
                        )

                    alternative_locations = []

                    if act["location"] in self.regions:
                        alternative_locations = [
                            act["location"],
                        ]

                    alternative_locations.extend(
                        [self.ecoinvent_to_iam_loc[act["location"]]]
                    )

                    for name_to_look_for, alt_loc in product(
                        names_to_look_for, alternative_locations
                    ):
                        if (
                            name_to_look_for,
                            exc[1],
                            alt_loc,
                            exc[-1],
                        ) in self.modified_datasets[
                            (self.model, self.scenario, self.year)
                        ][
                            "created"
                        ]:
                            entry = [(name_to_look_for, exc[1], alt_loc, exc[-1], 1.0)]

                            self.add_new_entry_to_cache(
                                location=act["location"],
                                exchange={
                                    "name": name_to_look_for,
                                    "product": exc[1],
                                    "location": exc[2],
                                    "unit": exc[-1],
                                },
                                allocated=[
                                    {
                                        "name": name_to_look_for,
                                        "product": exc[1],
                                        "location": alt_loc,
                                        "unit": exc[-1],
                                    }
                                ],
                                shares=[1.0],
                            )
                            break

                    if not entry:
                        if exc == (
                            act["name"],
                            act["reference product"],
                            act["location"],
                            act["unit"],
                        ):
                            entry = [
                                (
                                    act["name"],
                                    act["reference product"],
                                    act["location"],
                                    act["unit"],
                                    1.0,
                                )
                            ]

                    if not entry:
                        entry = [exc + (1.0,)]

                # summing up the amounts provided by the unwanted exchanges
                # and remove these unwanted exchanges from the dataset
                amount = sum(
                    e["amount"]
                    for e in excs_to_relink
                    if (e["name"], e["product"], e["location"], e["unit"]) == exc
                )

                if amount > 0:
                    new_exchanges.extend(
                        [
                            {
                                "name": e[0],
                                "product": e[1],
                                "amount": amount * e[-1],
                                "type": "technosphere",
                                "unit": e[3],
                                "location": e[2],
                            }
                            for e in entry
                        ]
                    )

            # make unique list of exchanges from new_exchanges
            # and sum the amounts of exchanges with the same name,
            # product, location and unit
            new_exchanges = [
                {
                    "name": name,
                    "product": prod,
                    "location": location,
                    "unit": unit,
                    "type": "technosphere",
                    "amount": sum([exc["amount"] for exc in exchanges]),
                }
                for (name, prod, location, unit), exchanges in groupby(
                    sorted(
                        new_exchanges,
                        key=itemgetter("name", "product", "location", "unit"),
                    ),
                    key=itemgetter("name", "product", "location", "unit"),
                )
            ]

            act["exchanges"] = [
                e
                for e in act["exchanges"]
                if (e["name"], e.get("product"), e.get("location"), e["unit"])
                not in [
                    (iex[0], iex[1], iex[2], iex[3]) for iex in unique_excs_to_relink
                ]
            ]
            act["exchanges"].extend(new_exchanges)

    def get_carbon_capture_rate(self, loc: str, sector: str) -> float:
        """
        Returns the carbon capture rate (between 0 and 1) as indicated by the IAM
        It is calculated as CO2 captured / (CO2 captured + CO2 emitted)

        :param loc: location of the dataset
        :param sector: name of the sector to look capture rate for
        :return: rate of carbon capture
        """

        if sector in self.iam_data.carbon_capture_rate.variables.values:
            rate = (
                self.iam_data.carbon_capture_rate.sel(
                    variables=sector,
                    region=loc,
                )
                .interp(year=self.year)
                .values
            )
        else:
            rate = 0

        return rate

    def create_ccs_dataset(
        self,
        loc: str,
        bio_co2_stored: float,
        bio_co2_leaked: float,
        sector: str = "cement",
    ) -> None:
        """
        Create a CCS dataset, reflecting the share of fossil vs. biogenic CO2.

        Source for CO2 capture and compression:
        https://www.sciencedirect.com/science/article/pii/S1750583613001230?via%3Dihub#fn0040

        :param loc: location of the dataset to create
        :param bio_co2_stored: share of biogenic CO2 over biogenic + fossil CO2
        :param bio_co2_leaked: share of biogenic CO2 leaked back into the atmosphere
        :param sector: name of the sector to look capture rate for
        :return: Does not return anything, but adds the dataset to the database.

        """

        # select the dataset
        # it is initially made for a cement plant, but it should be possible to
        # use it for any plant with a similar flue gas composition (CO2 concentration
        # and composition of the flue gas).
        dataset = ws.get_one(
            self.database,
            ws.equals(
                "name",
                "carbon dioxide, captured at cement production plant, "
                "with underground storage, post, 200 km",
            ),
            ws.equals("location", "RER"),
        )

        # duplicate the dataset
        ccs = wt.copy_to_new_location(dataset, loc)
        ccs["code"] = str(uuid.uuid4().hex)

        if sector != "cement":
            ccs["name"] = ccs["name"].replace("cement", sector)
            for e in ws.production(ccs):
                e["name"] = e["name"].replace("cement", sector)

        # relink the providers inside the dataset given the new location
        ccs = self.relink_technosphere_exchanges(ccs)

        if "input" in ccs:
            ccs.pop("input")

        # we first fix the biogenic CO2 permanent storage
        # this corresponds to the share of biogenic CO2
        # in the fossil + biogenic CO2 emissions of the plant

        for exc in ws.biosphere(
            ccs,
            ws.equals("name", "Carbon dioxide, in air"),
        ):
            exc["amount"] = bio_co2_stored

        if bio_co2_leaked > 0:
            # then the biogenic CO2 leaked during the capture process
            for exc in ws.biosphere(
                ccs,
                ws.equals("name", "Carbon dioxide, non-fossil"),
            ):
                exc["amount"] = bio_co2_leaked

        # the rest of CO2 leaked is fossil
        for exc in ws.biosphere(ccs, ws.equals("name", "Carbon dioxide, fossil")):
            exc["amount"] = 0.11 - bio_co2_leaked

        # we adjust the heat needs by subtraction 3.66 MJ with what
        # the plant is expected to produce as excess heat

        # Heat, as steam: 3.66 MJ/kg CO2 captured in 2020,
        # decreasing to 2.6 GJ/t by 2050, by looking at
        # the best-performing state-of-the-art technologies today
        # https://www.globalccsinstitute.com/wp-content/uploads/2022/05/State-of-the-Art-CCS-Technologies-2022.pdf
        # minus excess heat generated on site
        # the contribution of excess heat is assumed to be
        # 15% of the overall heat input with today's heat requirement
        # (see IEA 2018 cement roadmap report)

        heat_input = np.clip(np.interp(self.year, [2020, 2050], [3.66, 2.6]), 2.6, 3.66)
        excess_heat_generation = 3.66 * 0.15
        fossil_heat_input = heat_input - excess_heat_generation

        for exc in ws.technosphere(ccs, ws.contains("name", "steam production")):
            exc["amount"] = fossil_heat_input

        # then, we need to find local suppliers of electricity, water, steam, etc.
        ccs = self.relink_technosphere_exchanges(ccs)

        # add it to the list of created datasets
        self.modified_datasets[(self.model, self.scenario, self.year)][
            "created"
        ].append((ccs["name"], ccs["reference product"], ccs["location"], ccs["unit"]))

        # finally, we add this new dataset to the database
        self.database.append(ccs)

    def find_iam_efficiency_change(
        self,
        data: xr.DataArray,
        variable: Union[str, list],
        location: str,
    ) -> float:
        """
        Return the relative change in efficiency for `variable` in `location`
        relative to 2020.
        :param variable: IAM variable name
        :param location: IAM region
        :return: relative efficiency change (e.g., 1.05)
        """

        scaling_factor = (
            data.sel(region=location, variables=variable)
            .interp(year=self.year)
            .values.item(0)
        )

        if scaling_factor in (np.nan, np.inf):
            scaling_factor = 1

        return scaling_factor

    def write_log(self, dataset, status="created"):
        """
        Write log file.
        """

        logger.info(
            f"{status}|{self.model}|{self.scenario}|{self.year}|"
            f"{dataset['name']}|{dataset['location']}|"
        )

    def add_new_entry_to_cache(
        self,
        location: str,
        exchange: dict,
        allocated: List[dict],
        shares: List[float],
    ) -> None:
        """
        Add an entry to the cache.
        :param location: location
        :param exchange: exchange
        :param allocated: allocated
        :param shares: shares
        :return: cache
        """
        if "reference product" in exchange and "product" not in exchange:
            exchange["product"] = exchange["reference product"]

        exc_key = tuple(exchange[k] for k in ("name", "product", "location", "unit"))
        entry = [
            (
                e["name"],
                e["product"] if "product" in e else e["reference product"],
                e["location"],
                e["unit"],
                s,
            )
            for e, s in zip(allocated, shares)
        ]
        self.cache.setdefault(location, {}).setdefault(self.model, {})[exc_key] = entry

    def relink_technosphere_exchanges(
        self,
        dataset,
        exclusive=True,
        biggest_first=False,
        contained=True,
    ) -> dict:
        """Find new technosphere providers based on the location of the dataset.
        Designed to be used when the dataset's location changes, or when new datasets are added.
        Uses the name, reference product, and unit of the exchange to filter possible inputs. These must match exactly. Searches in the list of datasets ``data``.
        Will only search for providers contained within the location of ``dataset``, unless ``contained`` is set to ``False``, all providers whose location intersects the location of ``dataset`` will be used.
        A ``RoW`` provider will be added if there is a single topological face in the location of ``dataset`` which isn't covered by the location of any providing activity.
        If no providers can be found, `relink_technosphere_exchanes` will try to add a `RoW` or `GLO` providers, in that order, if available. If there are still no valid providers, a ``InvalidLink`` exception is raised, unless ``drop_invalid`` is ``True``, in which case the exchange will be deleted.
        Allocation between providers is done using ``allocate_inputs``; results seem strange if ``contained=False``, as production volumes for large regions would be used as allocation factors.
        Input arguments:
            * ``dataset``: The dataset whose technosphere exchanges will be modified.
            * ``data``: The list of datasets to search for technosphere product providers.
            * ``model``: The IAM model
            * ``exclusive``: Bool, default is ``True``. Don't allow overlapping locations in input providers.
            * ``drop_invalid``: Bool, default is ``False``. Delete exchanges for which no valid provider is available.
            * ``biggest_first``: Bool, default is ``False``. Determines search order when selecting provider locations. Only relevant is ``exclusive`` is ``True``.
            * ``contained``: Bool, default is ``True``. If true, only use providers whose location is completely within the ``dataset`` location; otherwise use all intersecting locations.
            * ``iam_regions``: List, lists IAM regions, if additional ones need to be defined.
        Modifies the dataset in place; returns the modified dataset."""

        new_exchanges = []
        technosphere = lambda x: x["type"] == "technosphere"

        list_loc = [k if isinstance(k, str) else k[1] for k in self.geo.geo.keys()] + [
            "RoW"
        ]

        for exc in filter(technosphere, dataset["exchanges"]):
            if (
                exc["name"],
                exc["product"],
                exc["location"],
                exc["unit"],
            ) in self.cache.get(dataset["location"], {}).get(self.model, {}):
                exchanges = self.cache[dataset["location"]][self.model][
                    (exc["name"], exc["product"], exc["location"], exc["unit"])
                ]

                if isinstance(exchanges, tuple):
                    exchanges = [exchanges]

                new_exchanges.extend(
                    [
                        {
                            "name": i[0],
                            "product": i[1],
                            "unit": i[3],
                            "location": i[2],
                            "type": "technosphere",
                            "amount": exc["amount"] * i[-1],
                        }
                        for i in exchanges
                    ]
                )

            else:
                kept = None
                key = (exc["name"], exc["product"], exc["unit"])
                possible_datasets = [
                    x for x in self.get_possibles(key) if x["location"] in list_loc
                ]

                possible_locations = [obj["location"] for obj in possible_datasets]

                _ = lambda x: (
                    x["name"],
                    x["reference product"]
                    if "reference product" in x
                    else x["product"],
                    x["location"],
                    x["unit"],
                )

                if len(possible_datasets) == 1:
                    assert (
                        possible_datasets[0]["name"] == exc["name"]
                    ), f"candidate: {_(possible_datasets[0])}, exc: {_(exc)}"
                    assert (
                        possible_datasets[0]["reference product"] == exc["product"]
                    ), f"candidate: {_(possible_datasets[0])}, exc: {_(exc)}"

                    self.cache.update(
                        {
                            dataset["location"]: {
                                self.model: {
                                    (
                                        exc["name"],
                                        exc["product"],
                                        exc["location"],
                                        exc["unit"],
                                    ): [
                                        (
                                            e["name"],
                                            e["product"],
                                            e["location"],
                                            e["unit"],
                                            s,
                                        )
                                        for e, s in zip(
                                            [new_exchange(exc, exc["location"], 1.0)],
                                            [1.0],
                                        )
                                    ]
                                }
                            }
                        }
                    )

                    new_exchanges.append(exc)
                    continue

                if dataset["location"] in possible_locations:
                    self.cache.update(
                        {
                            dataset["location"]: {
                                self.model: {
                                    (
                                        exc["name"],
                                        exc["product"],
                                        exc["location"],
                                        exc["unit"],
                                    ): [
                                        (
                                            e["name"],
                                            e["product"],
                                            e["location"],
                                            e["unit"],
                                            s,
                                        )
                                        for e, s in zip(
                                            [
                                                new_exchange(
                                                    exc, dataset["location"], 1.0
                                                )
                                            ],
                                            [1.0],
                                        )
                                    ]
                                }
                            }
                        }
                    )

                    exc["location"] = dataset["location"]
                    new_exchanges.append(exc)
                    continue

                if dataset["location"] in self.geo.iam_regions:
                    if any(
                        iloc in possible_locations
                        for iloc in self.geo.iam_to_ecoinvent_location(
                            dataset["location"]
                        )
                    ):
                        locs = [
                            iloc
                            for iloc in self.geo.iam_to_ecoinvent_location(
                                dataset["location"]
                            )
                            if iloc in possible_locations
                        ]

                        kept = [
                            ds for ds in possible_datasets if ds["location"] in locs
                        ]

                        if dataset["location"] == "World" and "GLO" in locs:
                            kept = [ds for ds in kept if ds["location"] == "GLO"]

                        allocated, share = allocate_inputs(exc, kept)
                        new_exchanges.extend(allocated)

                        self.add_new_entry_to_cache(
                            dataset["location"], exc, allocated, share
                        )
                        continue

                if not kept and any(
                    loc in possible_locations for loc in ["GLO", "RoW"]
                ):
                    kept = [
                        ds
                        for ds in possible_datasets
                        if ds["location"] in ["GLO", "RoW"]
                    ]
                    allocated, share = allocate_inputs(exc, kept)
                    new_exchanges.extend(allocated)
                    continue

                possible_locations = [
                    (self.model.upper(), p) if p in self.geo.iam_regions else p
                    for p in possible_locations
                ]

                if len(possible_datasets) > 0:
                    location = (
                        dataset["location"]
                        if dataset["location"] not in self.geo.iam_regions
                        else (self.model.upper(), dataset["location"])
                    )

                    gis_match = get_gis_match(
                        dataset,
                        location,
                        possible_locations,
                        self.geo,
                        contained,
                        exclusive,
                        biggest_first,
                    )

                    kept = [
                        ds
                        for loc in gis_match
                        for ds in possible_datasets
                        if ds["location"] == loc
                    ]

                    if kept:
                        missing_faces = self.geo.geo[location].difference(
                            set.union(*[self.geo.geo[obj["location"]] for obj in kept])
                        )
                        if missing_faces and "RoW" in possible_locations:
                            kept.extend(
                                [
                                    obj
                                    for obj in possible_datasets
                                    if obj["location"] == "RoW"
                                ]
                            )

                    if not kept and exc["name"].startswith("market group for"):
                        market_group_exc = exc.copy()
                        market_group_exc["name"] = market_group_exc["name"].replace(
                            "market group for", "market for"
                        )

                        key = (
                            market_group_exc["name"],
                            market_group_exc["product"],
                            market_group_exc["unit"],
                        )

                        possible_datasets = [
                            x
                            for x in self.get_possibles(key)
                            if x["location"] in list_loc
                        ]

                        possible_locations = [
                            obj["location"] for obj in possible_datasets
                        ]

                        possible_locations = [
                            (self.model.upper(), p) if p in self.geo.iam_regions else p
                            for p in possible_locations
                        ]

                        location = (
                            dataset["location"]
                            if dataset["location"] not in self.geo.iam_regions
                            else (self.model.upper(), dataset["location"])
                        )

                        gis_match = get_gis_match(
                            dataset,
                            location,
                            possible_locations,
                            self.geo,
                            contained,
                            exclusive,
                            biggest_first,
                        )

                        kept = [
                            ds
                            for loc in gis_match
                            for ds in possible_datasets
                            if ds["location"] == loc
                        ]

                        if kept:
                            exc = market_group_exc

                    kept = possible_datasets

                    allocated, share = allocate_inputs(exc, kept)
                    new_exchanges.extend(allocated)

                    # add to cache
                    self.add_new_entry_to_cache(
                        dataset["location"], exc, allocated, share
                    )

                else:
                    # there's no better candidate than the initial one
                    new_exchanges.append(exc)
                    # add to cache
                    self.add_new_entry_to_cache(
                        dataset["location"],
                        exc,
                        [exc],
                        [
                            1.0,
                        ],
                    )

        # make unique list of exchanges from new_exchanges
        # and sum the amounts of exchanges with the same name,
        # product, location and unit
        new_exchanges = [
            {
                "name": name,
                "product": prod,
                "location": location,
                "unit": unit,
                "type": "technosphere",
                "amount": sum([exc["amount"] for exc in exchanges]),
            }
            for (name, prod, location, unit), exchanges in groupby(
                sorted(
                    new_exchanges, key=itemgetter("name", "product", "location", "unit")
                ),
                key=itemgetter("name", "product", "location", "unit"),
            )
        ]

        dataset["exchanges"] = [
            exc for exc in dataset["exchanges"] if exc["type"] != "technosphere"
        ] + new_exchanges

        return dataset

    @lru_cache
    def get_possibles(self, key):
        """Filter a list of datasets ``data``,
        returning those with the save name,
        reference product, and unit as in ``exchange``.
        Returns a generator."""

        return [
            ds
            for ds in self.database
            if (ds["name"], ds["reference product"], ds["unit"]) == key
        ]
