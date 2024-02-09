"""
transformation.py contains the base class TransformationBase,
used by other classes (e.g. Transport, Electricity, Steel, Cement, etc.).
It provides basic methods usually used for electricity, cement, steel sector transformation
on the wurst database.
"""

import copy
import logging.config
import uuid
from collections import defaultdict
from collections.abc import ValuesView
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
from .filesystem_constants import DATA_DIR
from .geomap import Geomap
from .utils import get_fuel_properties

LOG_CONFIG = DATA_DIR / "utils" / "logging" / "logconfig.yaml"
# directory for log files
DIR_LOG_REPORT = Path.cwd() / "export" / "logs"
# if DIR_LOG_REPORT folder does not exist
# we create it
if not Path(DIR_LOG_REPORT).exists():
    Path(DIR_LOG_REPORT).mkdir(parents=True, exist_ok=True)

with open(LOG_CONFIG, encoding="utf-8") as f:
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
    exact_match: bool = False,
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

    if exact_match:
        filters = [
            ws.either(*[ws.equals("name", supplier) for supplier in names]),
        ]
    else:
        filters = [
            ws.either(*[ws.contains("name", supplier) for supplier in names]),
        ]

    filters += [
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

    def nonzero(x):
        return x if x != 0.0 else 1.0

    for dataset in dict_act:
        dict_act[dataset] /= nonzero(total_production_volume)

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

    def keep(x):
        return {
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
    return rescale_exchange(copied_exc, factor, remove_uncertainty=False)


def allocate_inputs(exc, lst):
    """
    Allocate the input exchanges in ``lst`` to ``exc``,
    using production volumes where possible, and equal splitting otherwise.
    Always uses equal splitting if ``RoW`` is present.
    """
    pvs = [reference_product(o).get("production volume", 0) for o in lst]

    if any((x > 0 for x in pvs)):
        # Allocate using production volume
        total = sum(pvs)
    else:
        # Allocate evenly
        total = len(lst)
        pvs = [1 for _ in range(total)]

    if lst[0]["name"] != exc["name"]:
        exc["name"] = lst[0]["name"]

    return (
        [
            new_exchange(exc, obj["location"], factor / total)
            for obj, factor in zip(lst, pvs)
            if factor > 0
        ],
        [p / total for p in pvs if p > 0],
    )


def filter_out_results(
    item_to_look_for: str, results: List[dict], field_to_look_at: str
) -> List[dict]:
    """Filters a list of results by a given field"""
    return [r for r in results if item_to_look_for not in r[field_to_look_at]]


def filter_technosphere_exchanges(exchanges: list):
    return filter(lambda x: x["type"] == "technosphere", exchanges)


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
        cache: dict = None,
        index: dict = None,
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
        self.heat_techs = mapping.generate_heat_map()
        self.system_model: str = system_model
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
        self.index = index or self.create_index()

    def create_index(self):
        idx = defaultdict(list)
        for ds in self.database:
            key = (copy.deepcopy(ds["name"]), copy.deepcopy(ds["reference product"]))
            idx[key].append(ds)
        return idx

    def add_to_index(self, ds: [dict, list, ValuesView]):
        if isinstance(ds, ValuesView):
            ds = list(ds)

        if isinstance(ds, dict):
            ds = [ds]

        for d in ds:
            key = (copy.deepcopy(d["name"]), copy.deepcopy(d["reference product"]))
            self.index[key].append(d)

    def remove_from_index(self, ds):
        key = (copy.deepcopy(ds["name"]), copy.deepcopy(ds["reference product"]))
        available_locations = [k["location"] for k in self.index[key]]
        if ds["location"] in available_locations:
            ds_to_remove = [
                d for d in self.index[key] if d["location"] == ds["location"]
            ][0]
            self.index[key].remove(ds_to_remove)

    def is_in_index(self, ds, location=None):
        if not any(key in ds for key in ["reference product", "product"]):
            raise KeyError(
                f"Dataset {ds['name']} does not have neither 'reference product' nor 'product' keys."
            )
        if "reference product" in ds:
            key = (ds["name"], ds["reference product"])
        else:
            key = (ds["name"], ds["product"])

        if location is None:
            return ds["location"] in [k["location"] for k in self.index[key]]

        return location in [k["location"] for k in self.index[key]]

    @lru_cache
    def select_multiple_suppliers(
        self,
        possible_names: Tuple[str],
        dataset_location: str,
        look_for: Tuple[str] = None,
        blacklist: Tuple[str] = None,
        exclude_region: Tuple[str] = None,
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
            "GLO",
            "Europe without Switzerland",
            "RER",
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

        if exclude_region:
            extra_filters.append(
                ws.exclude(
                    ws.either(*[ws.contains("location", x) for x in exclude_region])
                )
            )

        try:
            while not suppliers:
                suppliers = list(
                    ws.get_many(
                        self.database,
                        ws.either(
                            *[ws.contains("name", sup) for sup in possible_names]
                        ),
                        (
                            ws.either(
                                *[
                                    ws.equals("location", item)
                                    for item in possible_locations[counter]
                                ]
                            )
                            if isinstance(possible_locations[counter], list)
                            else ws.equals("location", possible_locations[counter])
                        ),
                        *extra_filters,
                    )
                )
                counter += 1
        except IndexError as err:

            suppliers = list(
                ws.get_many(
                    self.database,
                    ws.either(*[ws.contains("name", sup) for sup in possible_names]),
                    *extra_filters,
                )
            )

            if not suppliers:
                raise IndexError(
                    f"No supplier found for {possible_names} in {possible_locations}, "
                    f"looking for terms: {look_for} "
                    f"and with blacklist: {blacklist}"
                ) from err

        suppliers = get_shares_from_production_volume(suppliers)

        return suppliers

    def get_ecoinvent_locs(self) -> List[str]:
        """
        Rerun a list of unique locations in ecoinvent

        :return: list of locations
        :rtype: list
        """

        locs = list(set(a["location"] for a in self.database))

        # add Laos
        if "LA" not in locs:
            locs.append("LA")

        # add Fiji
        if "FJ" not in locs:
            locs.append("FJ")

        # add Guinea
        if "GN" not in locs:
            locs.append("GN")

        # add Guyana
        if "GY" not in locs:
            locs.append("GY")

        # add Sierra Leone
        if "SL" not in locs:
            locs.append("SL")

        # add Solomon Islands
        if "SB" not in locs:
            locs.append("SB")

        # add Uganda
        if "UG" not in locs:
            locs.append("UG")

        # add Afghanistan
        if "AF" not in locs:
            locs.append("AF")

        return locs

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
        if fuel_unit in ["kilogram", "cubic meter"]:
            try:
                lhv = self.fuels_specs[self.fuel_map_reverse[fuel_name]]["lhv"]

            except KeyError:
                lhv = 0
        elif fuel_unit == "kilowatt hour":
            lhv = 3.6
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

        # not_allowed = ["thermal"]
        # key = []
        # if "parameters" in dataset:
        #     key = list(
        #         key
        #         for key in dataset["parameters"]
        #         if "efficiency" in key and not any(item in key for item in not_allowed)
        #     )
        #
        # if len(key) > 0:
        #     if "parameters" in dataset:
        #         dataset["parameters"]["efficiency"] = dataset["parameters"][key[0]]
        #     else:
        #         dataset["parameters"] = {"efficiency": dataset["parameters"][key[0]]}
        #
        #     return dataset["parameters"][key[0]]

        energy_input = np.sum(
            np.sum(
                np.asarray(
                    [
                        self.calculate_input_energy(
                            exc["name"], exc["amount"], exc["unit"]
                        )
                        for exc in dataset["exchanges"]
                        if exc["name"] in fuel_filters
                        and exc["type"] == "technosphere"
                        and exc["amount"] > 0
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

        fallback_loc = None
        if "RoW" in d_map.values():
            fallback_loc = "RoW"
        else:
            if "GLO" in d_map.values():
                fallback_loc = "GLO"
            else:
                try:
                    fallback_loc = list(d_map.values())[0]
                except IndexError:
                    print(name, ref_prod, regions, d_map)

        return {region: d_map.get(region, fallback_loc) for region in regions}

    def fetch_proxies(
        self,
        name,
        ref_prod,
        production_variable=None,
        relink=True,
        regions=None,
        geo_mapping: dict = None,
        delete_original_dataset=False,
        empty_original_activity=True,
        exact_name_match=True,
        exact_product_match=False,
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
        :param delete_original_dataset: if True, delete original datasets from the database.
        :param empty_original_activity: if True, empty original activities from exchanges.
        :return: dictionary with IAM regions as keys, proxy datasets as values.
        """

        d_iam_to_eco = geo_mapping or self.region_to_proxy_dataset_mapping(
            name=name, ref_prod=ref_prod, regions=regions
        )

        d_act = {}

        ds_name, ds_ref_prod = [None, None]

        for region in d_iam_to_eco:
            # build filters
            if exact_name_match is True:
                filters = [
                    ws.equals("name", name),
                ]
            else:
                filters = [
                    ws.contains("name", name),
                ]
            if exact_product_match is True:
                filters.append(ws.equals("reference product", ref_prod))
            else:
                filters.append(ws.contains("reference product", ref_prod))

            filters.append(ws.equals("location", d_iam_to_eco[region]))

            try:
                dataset = ws.get_one(
                    self.database,
                    *filters,
                )
            except ws.MultipleResults as err:
                results = ws.get_many(
                    self.database,
                    *filters,
                )
                raise ws.MultipleResults(
                    err,
                    "A single dataset was expected, "
                    f"but found more than one for: {name, ref_prod}, : {[(r['name'], r['reference product'], r['location']) for r in results]}",
                )

            if not self.is_in_index(dataset, region):
                d_act[region] = copy.deepcopy(dataset)
                d_act[region]["location"] = region
                d_act[region]["code"] = str(uuid.uuid4().hex)

                for exc in ws.production(d_act[region]):
                    if "input" in exc:
                        del exc["input"]
                    if "location" in exc:
                        exc["location"] = region

                if "input" in d_act[region]:
                    del d_act[region]["input"]

                if production_variable is not None:
                    # Add `production volume` field
                    if isinstance(production_variable, str):
                        production_variable = [
                            production_variable,
                        ]

                    if isinstance(production_variable, list):
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

                    elif isinstance(production_variable, dict):
                        prod_vol = production_variable[region]
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
        for ds in ws.get_many(
            self.database,
            ws.equals("name", ds_name),
            ws.equals("reference product", ds_ref_prod),
        ):
            self.remove_from_index(ds)

        # empty original datasets
        # and make them link to new regional datasets
        if empty_original_activity is True:
            self.empty_original_datasets(
                name=ds_name,
                ref_prod=ds_ref_prod,
                loc_map=d_iam_to_eco,
                production_variable=production_variable,
                regions=regions,
            )

        if delete_original_dataset is True:
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
        production_variable: [str, dict],
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
        :param regions: regions to empty original datasets for
        :return: Does not return anything. Just empties the original dataset.
        """

        regions = regions or self.regions

        mapping = {}
        if loc_map:
            mapping = defaultdict(set)
            for v in loc_map.values():
                if self.geo.ecoinvent_to_iam_location(v) in loc_map.keys():
                    mapping[v].add(self.geo.ecoinvent_to_iam_location(v))
        existing_datasets = ws.get_many(
            self.database,
            ws.equals("name", name),
            ws.equals("reference product", ref_prod),
            ws.exclude(
                ws.either(*[ws.equals("location", loc) for loc in self.regions])
            ),
        )

        for existing_ds in existing_datasets:
            if existing_ds["location"] in mapping:
                locations = list(mapping[existing_ds["location"]])
            else:
                locations = [self.ecoinvent_to_iam_loc[existing_ds["location"]]]
                locations = [loc for loc in locations if loc in regions]

            if locations == [
                "World",
            ]:
                locations = [r for r in regions if r != "World"]

            if len(locations) == 0:
                continue

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

            def _(x):
                return x if x != 0.0 else 1.0

            if len(locations) == 1:
                existing_ds["exchanges"].append(
                    {
                        "name": existing_ds["name"],
                        "product": existing_ds["reference product"],
                        "amount": 1.0,
                        "unit": existing_ds["unit"],
                        "uncertainty type": 0,
                        "location": locations[0],
                        "type": "technosphere",
                    }
                )

            elif isinstance(production_variable, list) and all(
                i in self.iam_data.production_volumes.variables.values.tolist()
                for i in production_variable
            ):
                for location in locations:
                    share = (
                        self.iam_data.production_volumes.sel(
                            region=location, variables=production_variable
                        )
                        .interp(year=self.year)
                        .sum(dim="variables")
                        .values.item(0)
                    ) / _(
                        self.iam_data.production_volumes.sel(
                            region=locations, variables=production_variable
                        )
                        .interp(year=self.year)
                        .sum(dim=["variables", "region"])
                        .values.item(0)
                    )

                    if share > 0:
                        existing_ds["exchanges"].append(
                            {
                                "name": existing_ds["name"],
                                "product": existing_ds["reference product"],
                                "amount": share,
                                "unit": existing_ds["unit"],
                                "uncertainty type": 0,
                                "location": location,
                                "type": "technosphere",
                            }
                        )

            elif isinstance(production_variable, dict):
                existing_ds["exchanges"].extend(
                    [
                        {
                            "name": existing_ds["name"],
                            "product": existing_ds["reference product"],
                            "amount": share,
                            "unit": existing_ds["unit"],
                            "uncertainty type": 0,
                            "location": loc,
                            "type": "technosphere",
                        }
                        for loc, share in production_variable.items()
                    ]
                )

            else:
                share = 1 / len(locations)

                if share > 0:
                    existing_ds["exchanges"].extend(
                        [
                            {
                                "name": existing_ds["name"],
                                "product": existing_ds["reference product"],
                                "amount": share,
                                "unit": existing_ds["unit"],
                                "uncertainty type": 0,
                                "location": location,
                                "type": "technosphere",
                            }
                            for location in locations
                        ]
                    )

            self.remove_from_index(existing_ds)

            # add log
            self.write_log(dataset=existing_ds, status="empty")

    def relink_datasets(self, excludes_datasets=None, alt_names=None):
        """
        For a given exchange name, product, and unit, change its location to an IAM location,
        to effectively link to the newly built market(s)/activity(ies).
        :param excludes_datasets: list of datasets to exclude from relinking
        :param alt_names: list of alternative names to use for relinking
        """

        # Simplify default arguments
        alt_names = alt_names or []
        excludes_datasets = excludes_datasets or []

        for act in ws.get_many(
            self.database, ws.doesnt_contain_any("name", excludes_datasets)
        ):
            # Filter out exchanges to relink
            excs_to_relink = [
                e
                for e in ws.technosphere(act)
                if (not self.is_in_index(e) and e["amount"] != 0)
            ]

            if len(excs_to_relink) == 0:
                continue

            # make a dictionary with the names and amounts
            # of the technosphere exchanges to relink
            # to compare with the new exchanges
            excs_to_relink_dict = defaultdict(float)
            for exc in excs_to_relink:
                excs_to_relink_dict[exc["product"]] += exc["amount"]

            # Create a set of unique exchanges to relink
            # turn this into a list of dictionaries
            unique_excs_to_relink = [
                {
                    "name": exc["name"],
                    "product": exc["product"],
                    "location": exc["location"],
                    "unit": exc["unit"],
                }
                for exc in excs_to_relink
            ]
            # remove duplicates items in the list
            unique_excs_to_relink = [
                dict(t) for t in {tuple(d.items()) for d in unique_excs_to_relink}
            ]

            # Process exchanges to relink
            new_exchanges = self.process_exchanges_to_relink(
                act, unique_excs_to_relink, alt_names
            )

            # Update act["exchanges"] by removing the exchanges to relink
            act["exchanges"] = [e for e in act["exchanges"] if e not in excs_to_relink]
            # Update act["exchanges"] by adding new exchanges
            act["exchanges"].extend(new_exchanges)

            new_exchanges_dict = defaultdict(float)
            for exc in new_exchanges:
                new_exchanges_dict[exc["product"]] += exc["amount"]

            # compare with the original exchanges
            # if the amount is different, add a log
            for key in excs_to_relink_dict:
                assert (
                    key in new_exchanges_dict
                ), f"{key} not in {new_exchanges_dict} in dataset {act['name']}, {act['location']}"
                assert np.isclose(
                    excs_to_relink_dict[key],
                    new_exchanges_dict[key],
                    rtol=0.001,
                ), (
                    f"{excs_to_relink_dict[key]} != {new_exchanges_dict[key]} in dataset {act['name']}, {act['location']}."
                    f" Exchanges to relink: {excs_to_relink_dict}, new exchanges: {new_exchanges_dict}"
                )

    def process_exchanges_to_relink(self, act, unique_excs_to_relink, alt_names):
        new_exchanges = []
        for exc in unique_excs_to_relink:
            entries, amount = self.find_new_exchange_entries(act, exc, alt_names)
            if amount != 0:
                new_exchanges.extend(self.create_new_exchanges(entries, amount))
        # Make exchanges unique and sum amounts for duplicates
        return self.summarize_exchanges(new_exchanges)

    def get_exchange_from_cache(self, exc, loc):
        key = (
            exc["name"],
            exc["product"],
            exc["location"],
            exc["unit"],
        )

        return self.cache.get(loc, {}).get(self.model, {}).get(key)

    def find_alternative_locations(self, act, exc, alt_names):
        """
        Find alternative locations for an exchange, trying "market for" and "market group for"
        only if the initial search is unsuccessful.

        :param act: The activity dictionary.
        :param exc: A tuple representing the exchange (name, product, location, unit).
        :param alt_names: A list of alternative names to use for relinking.
        :return: A list of new exchange entries or an empty list if none are found.
        """
        names_to_look_for = [exc["name"]] + alt_names

        def allocate_exchanges(lst):
            """
            Allocate the input exchanges in ``lst`` to ``exc``,
            using production volumes where possible, and equal splitting otherwise.
            Always uses equal splitting if ``RoW`` is present.
            """

            pvs = []
            for o in lst:
                try:
                    ds = ws.get_one(
                        self.database,
                        ws.equals("name", o[0]),
                        ws.equals("reference product", o[1]),
                        ws.equals("location", o[2]),
                    )

                except ws.NoResults:
                    raise ws.NoResults(
                        f"Can't find {o[0]} {o[1]} {o[2]} in the database"
                    )

                for exc in ds["exchanges"]:
                    if exc["type"] == "production":
                        pvs.append(exc.get("production volume", 0))

            if any((x > 0 for x in pvs)):
                # Allocate using production volume
                total = sum(pvs)
            else:
                # Allocate evenly
                total = len(lst)
                pvs = [1 for _ in range(total)]

            return [p / total for p in pvs]

        # Function to search for new exchanges
        def search_for_new_exchanges(names):
            entries = []
            for name_to_look_for, alt_loc in product(
                set(names), set(alternative_locations)
            ):
                if (name_to_look_for, alt_loc) != (act["name"], act["location"]):
                    if self.is_in_index(
                        {
                            "name": name_to_look_for,
                            "product": exc["product"],
                            "location": alt_loc,
                            "unit": exc["unit"],
                        }
                    ):
                        entries.append(
                            (
                                name_to_look_for,
                                exc["product"],
                                alt_loc,
                                exc["unit"],
                                1.0,
                            )
                        )

            if len(entries) > 1 and any(
                x in ["World", "GLO", "RoW"] for x in [e[2] for e in entries]
            ):
                entries = [e for e in entries if e[2] not in ["World", "GLO", "RoW"]]

            if len(entries) > 1:
                shares = allocate_exchanges(entries)
                entries = [(e[0], e[1], e[2], e[3], s) for e, s in zip(entries, shares)]
                # remove entries that have a share of 0
                entries = [e for e in entries if e[-1] > 0]

            return entries

        # Start with the activity's location
        alternative_locations = [
            act["location"],
        ]

        # Add alternative locations based on mapping
        if act["location"] in self.ecoinvent_to_iam_loc:
            alternative_locations.append(self.ecoinvent_to_iam_loc[act["location"]])

        # Always include 'RoW', 'World', and 'GLO' as last resort options
        alternative_locations.extend(["RoW", "World", "GLO"])

        # Initial search with the provided names
        new_entries = search_for_new_exchanges(names_to_look_for)
        # check if the location of one of the entries matches with
        # the location of the activity to relink

        if len(new_entries) > 1:
            if any(e[2] == act["location"] for e in new_entries):
                new_entries = [e for e in new_entries if e[2] == act["location"]]
                # re-normalize the shares
                sum_shares = sum(e[-1] for e in new_entries)
                new_entries = [
                    (e[0], e[1], e[2], e[3], e[-1] / sum_shares) for e in new_entries
                ]

        if new_entries:
            return new_entries

        # If initial search fails, try with "market for" and "market group for"
        for prefix in ["market for", "market group for"]:
            if exc["name"].startswith(prefix):
                modified_name = exc["name"].replace(
                    prefix,
                    (
                        "market for"
                        if prefix == "market group for"
                        else "market group for"
                    ),
                )
                names_to_look_for.append(modified_name)

        # Second search with modified names
        return search_for_new_exchanges(names_to_look_for)

    def find_new_exchange_entries(self, act, exc, alt_names):
        entries = None

        if self.is_exchange_in_cache(exc, act["location"]):
            entries = self.get_exchange_from_cache(exc, act["location"])

        if not entries:
            entries = self.find_alternative_locations(act, exc, alt_names)

        if not entries:
            entries = [
                (exc["name"], exc["product"], exc["location"], exc["unit"]) + (1.0,)
            ]

        amount = sum(
            e["amount"]
            for e in ws.technosphere(act)
            if (e["name"], e["product"], e["location"], e["unit"])
            == (exc["name"], exc["product"], exc["location"], exc["unit"])
        )

        return entries, amount

    def create_new_exchanges(self, entries, amount):
        return [
            {
                "name": e[0],
                "product": e[1],
                "amount": amount * e[-1],
                "type": "technosphere",
                "unit": e[3],
                "location": e[2],
            }
            for e in entries
        ]

    def summarize_exchanges(self, new_exchanges):
        grouped_exchanges = groupby(
            sorted(
                new_exchanges, key=itemgetter("name", "product", "location", "unit")
            ),
            key=itemgetter("name", "product", "location", "unit"),
        )
        return [
            {
                "name": name,
                "product": prod,
                "location": loc,
                "unit": unit,
                "type": "technosphere",
                "amount": sum(e["amount"] for e in excs),
            }
            for (name, prod, loc, unit), excs in grouped_exchanges
        ]

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

        if not self.is_in_index(ccs):
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
            # 30% of heat requirement.

            heat_input = np.clip(
                np.interp(self.year, [2020, 2050], [3.66, 2.6]), 2.6, 3.66
            )
            excess_heat_generation = 0.3  # 30%
            fossil_heat_input = heat_input - (excess_heat_generation * heat_input)

            for exc in ws.technosphere(ccs, ws.contains("name", "steam production")):
                exc["amount"] = fossil_heat_input

            if sector != "cement":
                ccs["comment"] = ccs["comment"].replace("cement", sector)

            # then, we need to find local suppliers of electricity, water, steam, etc.
            ccs = self.relink_technosphere_exchanges(ccs)
            self.add_to_index(ccs)

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
        :param location: The location to which the cache entry corresponds.
        :param exchange: The exchange dictionary containing the data to cache.
        :param allocated: A list of dictionaries containing allocated exchanges.
        :param shares: A list of floats representing the shares for each allocated exchange.
        """
        # Ensure 'product' key is present in exchange.
        exchange.setdefault("product", exchange.get("reference product"))

        # Create a key for the cache entry.
        exc_key = (
            exchange["name"],
            exchange["product"],
            exchange["location"],
            exchange["unit"],
        )

        # Create the cache entry.
        entry = [
            (
                e.get("name", e.get("reference product")),
                e.get("product", e.get("reference product")),
                e["location"],
                e["unit"],
                s,
            )
            for e, s in zip(allocated, shares)
        ]

        # Initialize cache dictionary levels with setdefault.
        location_cache = self.cache.setdefault(location, {})
        model_cache = location_cache.setdefault(self.model, {})

        # Add the new entry to the cache.
        model_cache[exc_key] = entry

    def is_exchange_in_cache(self, exchange: dict, dataset_location: str) -> bool:
        """
        Check if an exchange is in the cache.
        :param exchange: The exchange dictionary to check.
        :param dataset_location: The location of the dataset.
        :return: True if the exchange is in the cache, False otherwise.

        """
        return (
            exchange["name"],
            exchange["product"],
            exchange["location"],
            exchange["unit"],
        ) in self.cache.get(dataset_location, {}).get(self.model, {})

    def process_cached_exchange(
        self, exchange: dict, dataset: dict, new_exchanges: list
    ) -> None:
        """
        Process a cached exchange. Adds the new exchanges to the list of new exchanges.
        :param exchange: The exchange dictionary to process.
        :param dataset: The dataset dictionary.
        :param new_exchanges: The list of new exchanges to add to the dataset.

        """
        exchanges = self.cache[dataset["location"]][self.model][
            (
                exchange["name"],
                exchange["product"],
                exchange["location"],
                exchange["unit"],
            )
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
                    "amount": exchange["amount"] * i[-1],
                }
                for i in exchanges
            ]
        )

    def process_uncached_exchange(
        self,
        exchange: dict,
        dataset: dict,
        new_exchanges: list,
        exclusive: bool,
        biggest_first: bool,
        contained: bool,
    ):
        """
        Process an uncached exchange. Adds the new exchanges to the list of new exchanges.
        :param exchange: The exchange dictionary to process.
        :param dataset: The dataset dictionary.
        :param new_exchanges: The list of new exchanges to add to the dataset.
        """

        # This function needs to handle the logic when
        # an exchange is not in the cache.
        key = (exchange["name"], exchange["product"])
        possible_datasets = self.index[key]

        if len(possible_datasets) == 0:
            if "market for" in exchange["name"]:
                key = (
                    exchange["name"].replace("market for", "market group for"),
                    exchange["product"],
                )
                possible_datasets = self.index[key]

        if len(possible_datasets) == 0:
            # search self.database for possible datasets
            possible_datasets = [
                ds
                for ds in self.database
                if ds["name"] == exchange["name"]
                and ds["reference product"] == exchange["product"]
            ]

        if len(possible_datasets) == 0:
            print(
                f"No possible datasets found for {key} in {dataset['name']} {dataset['location']}"
            )
            return [
                {
                    "name": exchange["name"],
                    "product": exchange["product"],
                    "unit": exchange["unit"],
                    "location": dataset["location"],
                    "type": "technosphere",
                    "amount": exchange["amount"],
                }
            ]

        if len(possible_datasets) == 1:
            self.handle_single_possible_dataset(
                exchange, possible_datasets, new_exchanges
            )

        else:
            self.handle_multiple_possible_datasets(
                exchange,
                dataset,
                possible_datasets,
                new_exchanges,
                exclusive,
                biggest_first,
                contained,
            )

    def handle_single_possible_dataset(
        self, exchange, possible_datasets, new_exchanges
    ):
        # If there's only one possible dataset, we can just use it
        single_dataset = possible_datasets[0]

        assert (
            single_dataset.get("reference product") == exchange["product"]
        ), f"Candidate: {single_dataset}, exchange: {exchange}"

        new_exc = exchange.copy()
        new_exc["location"] = single_dataset["location"]
        new_exc["name"] = single_dataset["name"]
        new_exc["product"] = single_dataset["reference product"]

        new_exchanges.append(new_exc)

    def new_exchange(self, exchange, location, amount_multiplier):
        # Create a new exchange dictionary with the modified location and amount
        return {
            "name": exchange["name"],
            "product": exchange["product"],
            "unit": exchange["unit"],
            "location": location,
            "type": "technosphere",
            "amount": exchange["amount"] * amount_multiplier,
        }

    def handle_multiple_possible_datasets(
        self,
        exchange: dict,
        dataset: dict,
        possible_datasets: list,
        new_exchanges: list,
        exclusive: bool,
        biggest_first: bool,
        contained: bool,
    ) -> None:
        # First, check if the dataset location itself is a possible match
        if dataset["location"] in [ds["location"] for ds in possible_datasets]:
            candidate = [
                ds for ds in possible_datasets if ds["location"] == dataset["location"]
            ][0]

            new_exc = exchange.copy()
            new_exc["location"] = candidate["location"]
            new_exc["name"] = candidate["name"]
            new_exc["product"] = candidate["reference product"]

            self.add_new_entry_to_cache(
                dataset["location"],
                exchange,
                [new_exc],
                [1.0],
            )

            new_exchanges.append(new_exc)
        else:
            # If more complex GIS matching or allocation is required,
            # we delegate to another function
            self.process_complex_matching_and_allocation(
                exchange,
                dataset,
                possible_datasets,
                new_exchanges,
                exclusive,
                biggest_first,
                contained,
            )

    def process_complex_matching_and_allocation(
        self,
        exchange: dict,
        dataset: dict,
        possible_datasets: list,
        new_exchanges: list,
        exclusive: bool,
        biggest_first: bool,
        contained: bool,
    ) -> None:
        # Check if the location of the dataset is within IAM regions
        if dataset["location"] in self.geo.iam_regions:
            self.handle_iam_region(exchange, dataset, possible_datasets, new_exchanges)

        elif dataset["location"] in ["GLO", "RoW", "World"]:
            # Handle global or rest-of-world scenarios
            self.handle_global_and_row_scenarios(
                exchange, dataset, possible_datasets, new_exchanges
            )

        else:
            # After the above checks, perform GIS matching if necessary
            self.perform_gis_matching(
                exchange,
                dataset,
                possible_datasets,
                new_exchanges,
                exclusive,
                biggest_first,
                contained,
            )

        # If there's still no match found, consider the default option
        self.handle_default_option(exchange, dataset, new_exchanges, possible_datasets)

    def handle_iam_region(self, exchange, dataset, possible_datasets, new_exchanges):
        # In IAM regions, we need to look for possible local datasets
        locs = [
            iloc
            for iloc in self.geo.iam_to_ecoinvent_location(dataset["location"])
            if iloc in [ds["location"] for ds in possible_datasets]
        ]

        if locs:
            kept = [ds for ds in possible_datasets if ds["location"] in locs]
            if dataset["location"] == "World" and "GLO" in locs:
                kept = [ds for ds in kept if ds["location"] == "GLO"]

            allocated, share = allocate_inputs(exchange, kept)

            new_exchanges.extend(allocated)
            self.add_new_entry_to_cache(dataset["location"], exchange, allocated, share)

    def handle_global_and_row_scenarios(
        self, exchange, dataset, possible_datasets, new_exchanges
    ):
        # Handle scenarios where the location is 'GLO' or 'RoW'
        possible_locations = [ds["location"] for ds in possible_datasets]
        if any(loc in possible_locations for loc in ["GLO", "RoW", "World"]):
            kept = [
                ds
                for ds in possible_datasets
                if ds["location"] in ["GLO", "RoW", "World"]
            ]
            allocated, share = allocate_inputs(exchange, kept)
            new_exchanges.extend(allocated)
            self.add_new_entry_to_cache(dataset["location"], exchange, allocated, share)

    def perform_gis_matching(
        self,
        exchange: dict,
        dataset: dict,
        possible_datasets: list,
        new_exchanges: list,
        exclusive: bool,
        biggest_first: bool,
        contained: bool,
    ) -> None:
        """
        Perform GIS matching for a dataset with a non-IAM location.

        :param exchange: The exchange dictionary to process.
        :param dataset: The dataset dictionary.
        :param possible_datasets: The list of possible datasets.
        :param new_exchanges: The list of new exchanges to add to the dataset.
        :param exclusive: Bool, default is ``True``. Don't allow overlapping locations in input providers.
        :param biggest_first: Bool, default is ``False``. Determines search order when selecting provider locations. Only relevant if ``exclusive`` is ``True``.
        :param contained: Bool, default is ``True``. If true, only use providers whose location is completely within the ``dataset`` location; otherwise use all intersecting locations.

        """
        # Perform GIS-based matching for location
        location = dataset["location"]
        # if regions contained in posisble location
        # we need to turn them into tuples (model, region)
        possible_locations = [ds["location"] for ds in possible_datasets]

        gis_match = self.get_gis_match(
            location,
            possible_locations,
            contained,
            exclusive,
            biggest_first,
        )

        kept = [
            ds for loc in gis_match for ds in possible_datasets if ds["location"] == loc
        ]

        if kept:
            allocated, share = allocate_inputs(exchange, kept)
            new_exchanges.extend(allocated)
            self.add_new_entry_to_cache(dataset["location"], exchange, allocated, share)

    def handle_default_option(
        self, exchange, dataset, new_exchanges, possible_datasets
    ):
        new_exc = None
        # Handle the default case where no better candidate is found
        if not self.is_exchange_in_cache(exchange, dataset["location"]):
            for default_location in ["RoW", "GLO", "World"]:
                if default_location in [x["location"] for x in possible_datasets]:
                    default_dataset = [
                        x
                        for x in possible_datasets
                        if x["location"] == default_location
                    ][0]

                    new_exc = exchange.copy()
                    new_exc["name"] = default_dataset["name"]
                    new_exc["product"] = default_dataset["reference product"]
                    new_exc["location"] = default_dataset["location"]
                    new_exchanges.append(new_exc)

                    break

        if new_exc is None and not self.is_exchange_in_cache(
            exchange, dataset["location"]
        ):
            new_exchanges.append(exchange)

    def find_candidates(
        self,
        dataset: dict,
        exclusive=True,
        biggest_first=False,
        contained=False,
    ):
        new_exchanges = []

        for exchange in filter_technosphere_exchanges(dataset["exchanges"]):
            if self.is_exchange_in_cache(exchange, dataset["location"]):
                self.process_cached_exchange(exchange, dataset, new_exchanges)
            else:
                self.process_uncached_exchange(
                    exchange,
                    dataset,
                    new_exchanges,
                    exclusive,
                    biggest_first,
                    contained,
                )

        return new_exchanges

    def relink_technosphere_exchanges(
        self,
        dataset,
        exclusive=True,
        biggest_first=False,
        contained=False,
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
            * ``biggest_first``: Bool, default is ``False``. Determines search order when selecting provider locations. Only relevant if ``exclusive`` is ``True``.
            * ``contained``: Bool, default is ``True``. If true, only use providers whose location is completely within the ``dataset`` location; otherwise use all intersecting locations.
            * ``iam_regions``: List, lists IAM regions, if additional ones need to be defined.
        Modifies the dataset in place; returns the modified dataset."""

        sum_before = sum(exc["amount"] for exc in dataset["exchanges"])

        # collect the name of exchange and the sum of amounts
        # as a dictionary, for all technosphere exchanges
        # in the dataset

        exchanges_before = defaultdict(float)
        for exc in dataset["exchanges"]:
            if exc["type"] == "technosphere":
                exchanges_before[exc["product"]] += exc["amount"]

        new_exchanges = self.find_candidates(
            dataset,
            exclusive=exclusive,
            biggest_first=biggest_first,
            contained=contained,
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
                "amount": sum(exc["amount"] for exc in exchanges),
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

        sum_after = sum(exc["amount"] for exc in dataset["exchanges"])

        assert np.allclose(sum_before, sum_after, rtol=1e-3), (
            f"Sum of exchanges before and after relinking is not the same: {sum_before} != {sum_after}"
            f"\n{dataset['name']}|{dataset['location']}"
        )

        # compare new exchanges with exchanges before
        exchanges_after = defaultdict(float)
        for exc in dataset["exchanges"]:
            if exc["type"] == "technosphere":
                exchanges_after[exc["product"]] += exc["amount"]

        assert set(exchanges_before.keys()) == set(exchanges_after.keys()), (
            f"Exchanges before and after relinking are not the same: {set(exchanges_before.keys())} != {set(exchanges_after.keys())}"
            f"\n{dataset['name']}|{dataset['location']}"
        )

        return dataset

    def get_gis_match(
        self,
        location,
        possible_locations,
        contained,
        exclusive,
        biggest_first,
    ):
        # prepare locations in possible_locations
        # all locations in possible_locations that are an IAM region
        # need to be converted to tuples with (model.upper(), location)
        # and other locations longer than 2 characters (other than GLO)
        # are converted to tuples with ("ecoinvent", location).

        possible_locations = [
            (
                (self.model.upper(), loc)
                if loc in self.regions
                else (
                    ("ecoinvent", loc)
                    if (len(loc) > 2 and loc not in ["GLO", "RoW"])
                    else loc
                )
            )
            for loc in possible_locations
        ]

        possible_locations = [loc for loc in possible_locations if loc in self.geo.geo]

        with resolved_row(possible_locations, self.geo.geo) as g:
            func = g.contained if contained else g.intersects

            gis_match = func(
                location,
                include_self=True,
                exclusive=exclusive,
                biggest_first=biggest_first,
                only=possible_locations,
            )

        return gis_match
