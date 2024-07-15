"""
Validates datapackages that contain external scenario data.
"""

import copy
import uuid
from collections import defaultdict
from typing import List, Union

import datapackage
import numpy as np
import pandas as pd
import wurst.searching as ws
import yaml
from datapackage import Package, exceptions, validate
from prettytable import PrettyTable
from schema import And, Optional, Schema, Use

from .geomap import Geomap
from .utils import load_constants

config = load_constants()


def find_iam_efficiency_change(
    variable: Union[str, list], location: str, efficiency_data, year: int
) -> float:
    """
    Return the relative change in efficiency for `variable` in `location`
    relative to 2020.
    :param variable: IAM variable name
    :param location: IAM region
    :return: relative efficiency change (e.g., 1.05)
    """

    scaling_factor = 1

    if variable in efficiency_data.variables.values:
        if year in efficiency_data.coords["year"].values:
            scaling_factor = (
                efficiency_data.sel(region=location, variables=variable, year=year)
            ).values.item(0)
        else:
            scaling_factor = (
                efficiency_data.sel(region=location, variables=variable).interp(
                    year=year
                )
            ).values.item(0)

        if scaling_factor in (np.nan, np.inf):
            scaling_factor = 1

    return scaling_factor


def flag_activities_to_adjust(
    dataset: dict,
    scenario_data: dict,
    year: int,
    dataset_vars: dict,
    region_proxy_mapping=None,
) -> dict:
    """
    Flag datasets that will need to be adjusted.
    :param dataset: dataset to be adjusted
    :param scenario_data: external scenario data
    :param year: year of the external scenario
    :param dataset_vars: variables of the dataset
    :return: dataset with additional info on variables to adjust
    """

    if "production volume variable" not in dataset_vars:
        regions = scenario_data["production volume"].region.values.tolist()
    else:
        data = scenario_data["production volume"].sel(
            variables=dataset_vars["production volume variable"]
        )
        # fetch regions which do not contain nan data
        regions = [
            r
            for r in data.region.values.tolist()
            if not np.isnan(data.sel(region=r).values).all()
        ]

    if "except regions" in dataset_vars:
        regions = [r for r in regions if r not in dataset_vars["except regions"]]

    dataset["regions"] = regions

    # add potential technosphere or biosphere filters
    if "efficiency" in dataset_vars:
        if len(dataset_vars["efficiency"]) > 0:
            dataset["adjust efficiency"] = True

            d_tech_filters = {
                k.get("variable"): [
                    k.get("includes").get("technosphere"),
                    {
                        region: find_iam_efficiency_change(
                            k["variable"],
                            region,
                            scenario_data["efficiency"],
                            year,
                        )
                        for region in regions
                    },
                ]
                for k in dataset_vars["efficiency"]
                if "technosphere" in k.get("includes", {})
            }

            d_tech_filters.update(
                {
                    k.get("variable"): [
                        None,
                        {
                            region: find_iam_efficiency_change(
                                k["variable"],
                                region,
                                scenario_data["efficiency"],
                                year,
                            )
                            for region in regions
                        },
                    ]
                    for k in dataset_vars["efficiency"]
                    if "includes" not in k
                }
            )

            d_bio_filters = {
                k.get("variable"): [
                    k.get("includes").get("biosphere"),
                    {
                        region: find_iam_efficiency_change(
                            k["variable"],
                            region,
                            scenario_data["efficiency"],
                            year,
                        )
                        for region in regions
                    },
                ]
                for k in dataset_vars["efficiency"]
                if "biosphere" in k.get("includes", {})
            }

            d_bio_filters.update(
                {
                    k.get("variable"): [
                        None,
                        {
                            region: find_iam_efficiency_change(
                                k["variable"],
                                region,
                                scenario_data["efficiency"],
                                year,
                            )
                            for region in regions
                        },
                    ]
                    for k in dataset_vars["efficiency"]
                    if "includes" not in k
                }
            )

            d_absolute_eff = {
                k.get("variable"): k.get("absolute", False)
                for k in dataset_vars["efficiency"]
            }

            if d_tech_filters:
                dataset["technosphere filters"] = d_tech_filters

            if d_bio_filters:
                dataset["biosphere filters"] = d_bio_filters

            if d_absolute_eff:
                dataset["absolute efficiency"] = d_absolute_eff

        # define exclusion filters
        for k in dataset_vars["efficiency"]:
            if "excludes" in k:
                if "technosphere" in k.get("excludes", {}):
                    dataset["excludes technosphere"] = {
                        k["variable"]: k["excludes"]["technosphere"]
                    }
                if "biosphere" in k.get("excludes", {}):
                    dataset["excludes biosphere"] = {
                        k["variable"]: k["excludes"]["biosphere"]
                    }

    if dataset_vars["replaces"]:
        dataset["replaces"] = dataset_vars["replaces"]

    if dataset_vars["replaces in"]:
        dataset["replaces in"] = dataset_vars["replaces in"]

    if dataset_vars["replacement ratio"] != 1.0:
        dataset["replacement ratio"] = dataset_vars["replacement ratio"]

    if dataset_vars["regionalize"]:
        dataset["regionalize"] = dataset_vars["regionalize"]

    if region_proxy_mapping is not None:
        dataset["region mapping"] = region_proxy_mapping

    if "production volume variable" in dataset_vars:
        dataset["production volume variable"] = dataset_vars[
            "production volume variable"
        ]

    return dataset


def check_inventories(
    configuration: dict,
    inventory_data: list,
    scenario_data: dict,
    database: list,
    year: int,
    model: str,
) -> tuple:
    """
    Check that the inventory data is valid.
    :param configuration: config file
    :param inventory_data: inventory data to check
    :param scenario_data: external scenario data
    :param database: database
    :param year: scenario year
    :param model: IAM model
    """

    # Initialize a defaultdict to count occurrences of each dataset
    dataset_usage = defaultdict(list)

    # Iterate over the production pathways
    for variable, pathway in configuration["production pathways"].items():
        if pathway["ecoinvent alias"].get("new dataset", False):
            continue
        # Extract the relevant keys for the dataset
        dataset_key = (
            pathway["ecoinvent alias"]["name"],
            pathway["ecoinvent alias"]["reference product"],
        )

        # Increment the usage count
        dataset_usage[dataset_key].append(variable)

    if any(len(x) > 1 for x in dataset_usage.values()):
        d = {}
        rows = []
        for k, v in dataset_usage.items():
            if len(v) > 1:
                rows.append((k[0][:50], k[1][:50]))
                for _, val in enumerate(v[1:]):
                    d[val] = k

        # print a Prettytable
        print("The following datasets will be duplicated:")
        table = PrettyTable()
        # adjust width of columns
        table.field_names = ["Name", "Reference product"]
        table._max_width = {"Name": 50, "Reference product": 50}
        for row in rows:
            table.add_row(row)
        print(table)

        for k, v in d.items():
            configuration["production pathways"][k]["ecoinvent alias"][
                "duplicate"
            ] = True
            configuration["production pathways"][k]["ecoinvent alias"][
                "name"
            ] += f"_{k}"

    geo = Geomap(model=model)

    d_datasets = {
        (
            val["ecoinvent alias"]["name"],
            val["ecoinvent alias"]["reference product"],
        ): {
            "exists in original database": val["ecoinvent alias"].get(
                "exists in original database", True
            ),
            "new dataset": val["ecoinvent alias"].get("new dataset", False),
            "duplicate": val["ecoinvent alias"].get("duplicate", False),
            "regionalize": val["ecoinvent alias"].get("regionalize", False),
            "mask": val["ecoinvent alias"].get("mask", None),
            "except regions": val.get(
                "except regions",
                [
                    # "World",
                ],
            ),
            "efficiency": val.get("efficiency", []),
            "replaces": val.get("replaces", []),
            "replaces in": val.get("replaces in", []),
            "replacement ratio": val.get("replacement ratio", 1),
            "production volume variable": val.get("production volume", {}).get(
                "variable"
            ),
        }
        for val in configuration["production pathways"].values()
    }

    # direct regionalization
    if "regionalize" in configuration:
        d_datasets.update(
            {
                (val["name"], val["reference product"]): {
                    "exists in original database": val.get(
                        "exists in original database", False
                    ),
                    "regionalize": True,
                    "new dataset": False,
                    "except regions": configuration["regionalize"].get(
                        "except regions", []
                    ),
                    "efficiency": val.get("efficiency", []),
                    "replaces": val.get("replaces", []),
                    "replaces in": val.get("replaces in", []),
                    "replacement ratio": val.get("replacement ratio", 1),
                }
                for val in configuration["regionalize"]["datasets"]
            }
        )

    list_datasets = [(i["name"], i["reference product"]) for i in inventory_data]

    try:
        assert all(
            (i[0], i[1]) in list_datasets
            for i, v in d_datasets.items()
            if not v["exists in original database"]
            and not v.get("new dataset")
            and not v.get("duplicate")
        )
    except AssertionError as e:
        list_missing_datasets = [
            i[0]
            for i, v in d_datasets.items()
            if not v["exists in original database"]
            and not v.get("new dataset")
            and (i[0], i[1]) not in list_datasets
        ]

        raise AssertionError(
            f"The following datasets are not in the inventory data: {list_missing_datasets}"
            f"\n Available datasets are: \n"
            f"{[list_datasets]}"
        ) from e

    # flag imported inventories
    for i, dataset in enumerate(inventory_data):
        key = (dataset["name"], dataset["reference product"])
        if key in d_datasets:
            if d_datasets[key]["exists in original database"] is False:
                dataset["custom scenario dataset"] = True
                data_vars = d_datasets[(dataset["name"], dataset["reference product"])]
                inventory_data[i] = flag_activities_to_adjust(
                    dataset, scenario_data, year, data_vars
                )

    def find_candidates_by_key(data, key):
        """Filter data for items matching the key (name and reference product)."""
        return [
            item
            for item in data
            if item["name"] == key[0] and item["reference product"] == key[1]
        ]

    def filter_candidates_by_mask(candidates, mask):
        """Exclude candidates containing the mask in their name."""
        return [candidate for candidate in candidates if mask not in candidate["name"]]

    def identify_potential_candidates(database, inventory_data, key, mask):
        """Identify and return potential candidates based on key and mask."""
        candidates = find_candidates_by_key(database + inventory_data, key)
        if mask:
            candidates = filter_candidates_by_mask(candidates, mask)
        return candidates

    def adjust_candidate(candidate, scenario_data, year, val, region=None):
        """Adjust a single candidate with scenario data."""
        flag_activities_to_adjust(candidate, scenario_data, year, val, region)

    def handle_single_candidate(candidates, scenario_data, year, val):
        """Handle case where there is exactly one candidate."""
        adjust_candidate(candidates[0], scenario_data, year, val)

    def perform_region_checks(potential_candidates, scenario_data):
        """
        Perform geographic region checks to shortlist candidates.

        This includes checking for direct matches, containment, intersection,
        and mapping between different geographic naming conventions.
        """
        # Create a shortlist dictionary with None as initial values
        short_listed = {
            r: None for r in scenario_data["production volume"].region.values
        }

        # Convert potential_candidates to a dictionary for O(1) lookups by location
        candidates_by_location = {c["location"]: c for c in potential_candidates}

        # Sort candidate locations once, considering fallback_locations
        fallback_locations = ["GLO", "RoW"]
        sorted_candidate_locations = sorted(
            candidates_by_location.keys(),
            key=lambda x: x in fallback_locations,
        )

        reasons = {}

        # Function to check and assign the first matching candidate
        def assign_candidate_if_empty(region, loc):
            if short_listed[region] is None and loc in candidates_by_location:
                short_listed[region] = candidates_by_location[loc]

        # Define check functions
        def direct_match(region, location):
            if region == location:
                reasons[region] = "direct match"
            return region == location

        def iam_match(region, location):
            if location in geo.iam_regions:
                reasons[region] = "IAM match"
                return region in geo.iam_to_ecoinvent_location(location)
            return False

        def contained_match(region, location):
            if location not in geo.iam_regions:
                try:
                    if region in geo.geo.contained(location):
                        reasons[region] = "contained match"
                    return region in geo.geo.contained(location)
                except KeyError:
                    return False
            return False

        def intersects_match(region, location):
            if location not in geo.iam_regions:
                try:
                    if region in geo.geo.intersects(location):
                        reasons[region] = "intersects match"
                    return region in geo.geo.intersects(location)
                except KeyError:
                    return False
            return False

        def ecoinvent_match(region, location):
            if geo.ecoinvent_to_iam_location(location) == "World":
                return False

            if region == geo.ecoinvent_to_iam_location(location):
                reasons[region] = "ecoinvent to IAM match"
            return region == geo.ecoinvent_to_iam_location(location)

        def fallback_match(region, location):
            if location in fallback_locations:
                reasons[region] = "fallback match"
            return location in fallback_locations

        # Ordered list of check functions
        check_functions = [
            direct_match,
            iam_match,
            contained_match,
            intersects_match,
            ecoinvent_match,
            fallback_match,
        ]

        # Perform checks in order of preference
        for region in short_listed:
            for location in sorted_candidate_locations:
                for check_func in check_functions:
                    if check_func(region, location):
                        # check the dataset was not previously emptied
                        if (
                            candidates_by_location[location].get("emptied", False)
                            is False
                        ):
                            assign_candidate_if_empty(region, location)
                            break

            if short_listed[region] is None:
                if "RoW" in sorted_candidate_locations:
                    assign_candidate_if_empty(region, "RoW")
                else:
                    assign_candidate_if_empty(region, location)

        return short_listed

    def short_list_candidates(candidates, scenario_data):
        """Shortlist candidates based on region logic."""
        short_listed = perform_region_checks(candidates, scenario_data)
        return short_listed

    def adjust_candidates_or_raise_error(
        candidates,
        scenario_data,
        key,
        year,
        val,
        inventory_data,
    ) -> [None, list]:
        """Adjust candidates if possible or raise an error if no valid candidates are found."""
        if not candidates:
            if not find_candidates_by_key(inventory_data, key):
                raise ValueError(
                    f"Dataset {key[0]} and {key[1]} is not found in the original database."
                )
            return  # Skip further processing if no candidates found.

        if len(candidates) == 1:
            handle_single_candidate(candidates, scenario_data, year, val)
            return candidates
        else:
            short_listed = short_list_candidates(candidates, scenario_data)
            for region, ds in short_listed.items():
                if ds is not None:
                    adjust_candidate(
                        ds,
                        scenario_data,
                        year,
                        val,
                        {r: d["location"] for r, d in short_listed.items()},
                    )
                else:
                    print(f"No candidate found for {key[0]} and {key[1]} for {region}.")
            return list(short_listed.values())

    for key, val in d_datasets.items():
        if val.get("exists in original database"):
            mask = val.get("mask")

            duplicate_name = None
            if val.get("duplicate") is True:
                duplicate_name = key[0]
                key = (key[0].split("_")[0], key[1])

            potential_candidates = identify_potential_candidates(
                database, inventory_data, key, mask
            )
            candidates = adjust_candidates_or_raise_error(
                potential_candidates,
                scenario_data,
                key,
                year,
                val,
                inventory_data,
            )

            if duplicate_name:
                for candidate in candidates:
                    # deep copy the candidate
                    ds = copy.deepcopy(candidate)
                    ds["code"] = str(uuid.uuid4().hex)
                    ds["name"] = duplicate_name
                    for exc in ws.production(ds):
                        exc["name"] = duplicate_name
                    database.append(ds)

    return inventory_data, database, configuration


def check_datapackage(datapackage: datapackage.Package):
    # validate package descriptor

    try:
        validate(datapackage.descriptor)
    except exceptions.ValidationError as exception:
        raise exception

    if "config" in [i.name for i in datapackage.resources] and "scenario_data" not in [
        i.name for i in datapackage.resources
    ]:
        raise ValueError(
            "If the resource 'config' is present in the datapackage,"
            "so must the resource 'scenario_data'."
        )

    if "scenario_data" in [i.name for i in datapackage.resources] and "config" not in [
        i.name for i in datapackage.resources
    ]:
        raise ValueError(
            "If the resource 'scenario_data' is present in the datapackage,"
            " so must the resource 'config'."
        )

    assert (
        datapackage.descriptor["ecoinvent"]["version"]
        in config["SUPPORTED_EI_VERSIONS"]
    ), f"The ecoinvent version in the datapackage is not supported. Must be one of {config['SUPPORTED_EI_VERSIONS']}."

    if (
        sum(
            s.name == y.name
            for s in datapackage.resources
            for y in datapackage.resources
        )
        / len(datapackage.resources)
        > 1
    ):
        raise ValueError(f"Two or more resources in the datapackage are similar.")


def list_all_iam_regions(configuration):
    """
    List all IAM regions in the config file.
    :param configuration: config file
    :return: list of IAM regions
    """

    list_regions = []

    for k, v in configuration.items():
        if k.startswith("LIST_"):
            list_regions.extend(v)

    return list_regions


def check_config_file(datapackage: datapackage.Package) -> int:

    resource = datapackage.get_resource("config")
    config_file = yaml.safe_load(resource.raw_read())

    file_schema = Schema(
        {
            "production pathways": {
                str: {
                    "production volume": {
                        "variable": str,
                    },
                    "ecoinvent alias": {
                        "name": str,
                        "reference product": str,
                        Optional("exists in original database"): bool,
                        Optional("mask"): str,
                        Optional("new dataset"): bool,
                        Optional("regionalize"): bool,
                        Optional("ratio"): float,
                    },
                    Optional("efficiency"): [
                        {
                            "variable": str,
                            Optional("reference year"): And(
                                Use(int), lambda n: 2005 <= n <= 2100
                            ),
                            Optional("includes"): {
                                Optional("technosphere"): list,
                                Optional("biosphere"): list,
                            },
                            Optional("excludes"): {
                                Optional("technosphere"): list,
                                Optional("biosphere"): list,
                            },
                            Optional("absolute"): bool,
                        }
                    ],
                    Optional("except regions"): And(
                        list,
                        Use(list),
                        lambda s: all(i in list_all_iam_regions(config) for i in s),
                    ),
                    Optional("replaces"): [
                        {
                            "name": str,
                            "product": str,
                            Optional("location"): str,
                            Optional("operator"): str,
                        }
                    ],
                    Optional("replaces in"): list,
                    Optional("replacement ratio"): float,
                },
            },
            Optional("markets"): [
                {
                    "name": str,
                    "reference product": str,
                    "unit": str,
                    "includes": And(
                        list,
                        Use(list),
                        lambda s: all(
                            i in config_file["production pathways"] for i in s
                        ),
                    ),
                    Optional("add"): [
                        {
                            Optional("name"): str,
                            Optional("reference product"): str,
                            Optional("categories"): str,
                            Optional("unit"): str,
                            Optional("amount"): float,
                        }
                    ],
                    Optional("except regions"): And(
                        list,
                        Use(list),
                        lambda s: all(
                            i
                            in config["LIST_REMIND_REGIONS"]
                            + config["LIST_IMAGE_REGIONS"]
                            for i in s
                        ),
                    ),
                    Optional("replaces"): [
                        {
                            "name": str,
                            "product": str,
                            Optional("location"): str,
                            Optional("operator"): str,
                        }
                    ],
                    Optional("replaces in"): list,
                    Optional("is fuel"): dict,
                    Optional("replacement ratio"): float,
                    Optional("waste market"): bool,
                    Optional("efficiency"): [
                        {
                            "variable": str,
                            Optional("reference year"): And(
                                Use(int), lambda n: 2005 <= n <= 2100
                            ),
                            Optional("includes"): {
                                Optional("technosphere"): list,
                                Optional("biosphere"): list,
                            },
                            Optional("excludes"): {
                                Optional("technosphere"): list,
                                Optional("biosphere"): list,
                            },
                        }
                    ],
                }
            ],
            Optional("regionalize"): {
                "datasets": [
                    {
                        "name": str,
                        "reference product": str,
                        Optional("exists in original database"): bool,
                    }
                ],
                Optional("except regions"): And(
                    list,
                    Use(list),
                    lambda s: all(
                        i
                        in config["LIST_REMIND_REGIONS"] + config["LIST_IMAGE_REGIONS"]
                        for i in s
                    ),
                ),
            },
        }
    )

    file_schema.validate(config_file)

    if "markets" in config_file:
        # check that providers composing the market
        # are listed

        for market in config_file["markets"]:
            try:
                [
                    (
                        config_file["production pathways"][a]["ecoinvent alias"][
                            "name"
                        ],
                        config_file["production pathways"][a]["ecoinvent alias"][
                            "reference product"
                        ],
                    )
                    for a in market["includes"]
                ]
            except KeyError as err:
                raise ValueError(
                    "One of more providers listed under `markets/includes` is/are not listed "
                    "under `production pathways`."
                ) from err

    needs_imported_inventories = False

    resource = datapackage.get_resource("config")
    config_file = yaml.safe_load(resource.raw_read())

    if len(list(config_file["production pathways"].keys())) != sum(
        get_recursively(
            config_file["production pathways"], "exists in original database"
        )
    ):
        needs_imported_inventories = True

    return needs_imported_inventories


def check_scenario_data_file(
    datapackage: datapackage.Package, scenario: str
) -> Package:

    scenarios = datapackage.descriptor["scenarios"]
    resource = datapackage.get_resource("scenario_data")
    scenario_data = resource.read()
    scenario_headers = resource.headers

    df = pd.DataFrame(scenario_data, columns=scenario_headers)

    resource = datapackage.get_resource("config")
    config_file = yaml.safe_load(resource.raw_read())

    mandatory_fields = [
        "scenario",
        "region",
        "variables",
        "unit",
    ]
    if not all(v in df.columns for v in mandatory_fields):
        raise ValueError(
            f"One or several mandatory column are missing "
            f"in the scenario data file. Mandatory columns: {mandatory_fields}."
        )

    years_cols = []
    for header in scenario_headers:
        try:
            years_cols.append(int(header))
        except ValueError:
            continue

    if not all(2005 <= y <= 2100 for y in years_cols):
        raise ValueError(
            f"One or several of the years provided in the scenario data file are "
            "out of boundaries (2005 - 2100)."
        )

    if len(pd.isnull(df).sum()[pd.isnull(df).sum() > 0]) > 0:
        raise ValueError(
            f"The following columns in the scenario data "
            f"contains empty cells.\n{pd.isnull(df).sum()[pd.isnull(df).sum() > 0]}."
        )

    available_scenarios = df["scenario"].unique()

    if not all(
        s in available_scenarios for s in scenarios
    ):  # check that all scenarios are available in the scenario file
        print(
            "The following scenarios listed in the json file "
            "are not available in the scenario data file:"
        )
        print(set(s for s in scenarios if s not in available_scenarios))
        raise ValueError(
            f"One or several scenarios are not available in the scenario file."
        )

    if scenario not in available_scenarios:
        raise ValueError(
            f"The scenario {scenario} is not available in the scenario file."
            f"Available scenarios are: {available_scenarios}."
        )

    if not all(
        v in df["variables"].unique() for v in get_recursively(config_file, "variable")
    ):
        list_unfound_variables = [
            p
            for p in get_recursively(config_file, "variable")
            if p not in df["variables"].unique()
        ]

        raise ValueError(
            "The following variables from the configuration file "
            f"cannot be found in the scenario file.: {list_unfound_variables}"
        )

    if not all(
        v in df["variables"].unique() for v in get_recursively(config_file, "variable")
    ):
        missing_variables = [
            v
            for v in get_recursively(config_file, "variable")
            if v not in df["variables"].unique()
        ]
        raise ValueError(
            f"One or several variable names in the configuration file "
            f"cannot be found in the scenario data file: {missing_variables}."
        )

    try:
        np.array_equal(df.iloc[:, 6:], df.iloc[:, 6:].astype(float))
    except ValueError as e:
        raise TypeError(
            f"All values provided in the time series must be numerical "
            f"in the scenario data file."
        ) from e

    return datapackage


def get_recursively(search_dict, field):
    """Takes a dict with nested lists and dicts,
    and searches all dicts for a key of the field
    provided.
    """
    fields_found = []

    for key, value in search_dict.items():
        if key == field:
            fields_found.append(value)

        elif isinstance(value, dict):
            results = get_recursively(value, field)
            for result in results:
                fields_found.append(result)

        elif isinstance(value, list):
            for item in value:
                if isinstance(item, dict):
                    more_results = get_recursively(item, field)
                    for another_result in more_results:
                        fields_found.append(another_result)

    return fields_found


def check_external_scenarios(external_scenarios: list) -> list:
    """
    Check external scenarios.
    :param external_scenarios: external scenario data

    :return: external scenario data, IAM scenario data
    """

    for external_scenario in external_scenarios:
        # Validate datapackage
        check_datapackage(external_scenario["data"])

        # Validate yaml config file
        check_config_file(external_scenario["data"])

        # Validate scenario data
        external_scenario["data"] = check_scenario_data_file(
            external_scenario["data"], external_scenario["scenario"]
        )

    return external_scenarios


def fetch_dataset_description_from_production_pathways(
    configuration: dict, item: str
) -> tuple:
    for p, v in configuration["production pathways"].items():
        if p == item:
            if "exists in original database" not in v["ecoinvent alias"]:
                v["ecoinvent alias"].update({"exists in original database": True})

            if "new dataset" not in v["ecoinvent alias"]:
                v["ecoinvent alias"].update({"new dataset": False})

            return (
                v["ecoinvent alias"]["name"],
                v["ecoinvent alias"]["reference product"],
                v["ecoinvent alias"]["exists in original database"],
                v["ecoinvent alias"]["new dataset"],
            )
    return None, None, None, None
