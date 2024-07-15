"""
export.py contains all the functions to format, prepare and export databases.
"""

import csv
import datetime
import json
import multiprocessing as mp
import os
import re
import uuid
from collections import defaultdict
from datetime import datetime
from functools import lru_cache
from multiprocessing.pool import ThreadPool as Pool
from pathlib import Path
from typing import Any, Dict, List

import numpy as np
import pandas as pd
import sparse
import yaml
from datapackage import Package
from pandas import DataFrame
from prettytable import PrettyTable
from scipy import sparse as nsp

from . import __version__
from .data_collection import get_delimiter
from .filesystem_constants import DATA_DIR
from .inventory_imports import get_correspondence_bio_flows
from .utils import dump_database, load_database, reset_all_codes
from .validation import BaseDatasetValidator

FILEPATH_SIMAPRO_UNITS = DATA_DIR / "utils" / "export" / "simapro_units.yml"
FILEPATH_SIMAPRO_COMPARTMENTS = (
    DATA_DIR / "utils" / "export" / "simapro_compartments.yml"
)
CORRESPONDENCE_BIO_FLOWS = (
    DATA_DIR / "utils" / "export" / "correspondence_biosphere_flows.yaml"
)

# current working directory
DIR_DATAPACKAGE = Path.cwd() / "export" / "datapackage"
DIR_DATAPACKAGE_TEMP = Path.cwd() / "export" / "temp"


def get_simapro_units() -> Dict[str, str]:
    """
    Load a dictionary that maps brightway2 unit to Simapro units.
    :return: a dictionary that maps brightway2 unit to Simapro units
    """

    with open(FILEPATH_SIMAPRO_UNITS, "r", encoding="utf-8") as stream:
        simapro_units = yaml.safe_load(stream)

    return simapro_units


def get_simapro_compartments() -> Dict[str, str]:
    """
    Load a dictionary that maps brightway2 unit to Simapro compartments.
    :return: a dictionary that maps brightway2 unit to Simapro compartments.
    """

    with open(FILEPATH_SIMAPRO_COMPARTMENTS, "r", encoding="utf-8") as stream:
        simapro_comps = yaml.safe_load(stream)

    return simapro_comps


def load_simapro_categories():
    """Load a dictionary with categories to use for Simapro export"""

    # Load the matching dictionary
    filename = "simapro_classification.csv"
    filepath = DATA_DIR / "utils" / "export" / filename
    if not filepath.is_file():
        raise FileNotFoundError(
            "The dictionary of Simapro categories could not be found."
        )
    with open(filepath, encoding="latin1") as f:
        csv_list = [[val.strip() for val in r.split(";")] for r in f.readlines()]
    _, *data = csv_list

    dict_cat = {}
    for row in data:
        _, cat_code, category_1, category_2, category_3 = row
        dict_cat[str(cat_code)] = {
            "category 1": category_1,
            "category 2": category_2,
            "category 3": category_3,
        }

    return dict_cat


def get_simapro_category_of_exchange():
    """Load a dictionary with categories to use for Simapro export based on ei 3.7"""

    # Load the matching dictionary
    filename = "simapro_categories.csv"
    filepath = DATA_DIR / "utils" / "export" / filename
    if not filepath.is_file():
        raise FileNotFoundError(
            "The dictionary of Simapro categories could not be found."
        )
    with open(
        filepath,
        encoding="latin-1",
    ) as file:
        csv_reader = csv.DictReader(file)

        dict_cat = {}
        # first row are headers
        for row in csv_reader:
            dict_cat[(row["name"].lower(), row["product"].lower())] = dict(row)

    return dict_cat


def load_references():
    """Load a dictionary with references of datasets"""

    # Load the matching dictionary
    filename = "references.csv"
    filepath = DATA_DIR / "utils" / "export" / filename
    if not filepath.is_file():
        raise FileNotFoundError("The dictionary of references could not be found.")
    with open(filepath, encoding="utf-8") as file:
        csv_list = [[val.strip() for val in r.split(";")] for r in file.readlines()]
    _, *data = csv_list

    dict_reference = {}
    for row in data:
        name, source, description = row
        dict_reference[name] = {"source": source, "description": description}

    return dict_reference


def get_simapro_biosphere_dictionnary():
    """
    Load a dictionary with biosphere flows to use for Simapro export.
    """

    # Load the matching dictionary between ecoinvent and Simapro biosphere flows
    filename = "simapro-biosphere.json"
    filepath = DATA_DIR / "utils" / "export" / filename
    if not filepath.is_file():
        raise FileNotFoundError(
            "The dictionary of biosphere flow match between ecoinvent "
            "and Simapro could not be found."
        )
    with open(filepath, encoding="utf-8") as json_file:
        data = json.load(json_file)
    dict_bio = {}
    for row in data:
        dict_bio[row[2]] = row[1]

    return dict_bio


def create_index_of_A_matrix(database):
    """
    Create a dictionary with row/column indices of the exchanges
    matrix as key and a tuple (activity name, reference product,
    unit, location) as value.
    :return: a dictionary to map indices to activities
    :rtype: dict
    """
    return {
        (
            database[i]["name"],
            database[i]["reference product"],
            database[i]["unit"],
            database[i]["location"],
        ): i
        for i in range(0, len(database))
    }


def rev_index(inds: dict) -> dict:
    """
    Reverse the index of the A matrix.
    """
    return {v: k for k, v in inds.items()}


def create_codes_index_of_exchanges_matrix(database):
    """
    Create a dictionary with row/column indices of the A matrix
    as key and the activity code as value.
    :return: a dictionary to map indices to activity codes
    :rtype: dict
    """

    return {database[i]["code"]: i for i in range(0, len(database))}


def create_codes_index_of_biosphere_flows_matrix(version):
    """
    Create a dictionary with row/column indices of the biosphere matrix
    """
    data = biosphere_flows_dictionary(version)

    return {v: k for k, v in enumerate(data.values())}


def create_index_of_biosphere_flows_matrix(version):
    data = biosphere_flows_dictionary(version)

    return {v: k for k, v in enumerate(data.keys())}


def create_codes_and_names_of_tech_matrix(database: List[dict]):
    """
    Create a dictionary a tuple (activity name, reference product,
    unit, location) as key, and its code as value.
    :return: a dictionary to map indices to activities
    :rtype: dict
    """
    return {
        (
            i["name"],
            i["reference product"],
            i["unit"],
            i["location"],
        ): i["code"]
        for i in database
    }


@lru_cache
def biosphere_flows_dictionary(version):
    """
    Create a dictionary with biosphere flows
    (name, category, sub-category, unit) -> code
    """
    if version == "3.9":
        fp = DATA_DIR / "utils" / "export" / "flows_biosphere_39.csv"
    elif version == "3.10":
        fp = DATA_DIR / "utils" / "export" / "flows_biosphere_310.csv"
    else:
        fp = DATA_DIR / "utils" / "export" / "flows_biosphere_38.csv"

    if not Path(fp).is_file():
        raise FileNotFoundError("The dictionary of biosphere flows could not be found.")

    csv_dict = {}

    with open(fp, encoding="utf-8") as file:
        input_dict = csv.reader(
            file,
            delimiter=get_delimiter(filepath=fp),
        )
        for row in input_dict:
            csv_dict[(row[0], row[1], row[2], row[3])] = row[-1]

    return csv_dict


def get_list_unique_acts(scenarios: List[dict]) -> list:
    """
    Get a list of unique activities from a list of databases
    :param scenarios: list of databases
    :return: list of unique activities
    """

    list_unique_acts = []
    for db in scenarios:
        for ds in db["database"]:
            list_unique_acts.extend(
                [
                    (
                        a["name"],
                        a.get("product"),
                        a.get("categories"),
                        a.get("location"),
                        a["unit"],
                        a["type"],
                    )
                    for a in ds["exchanges"]
                ]
            )
    return list(set(list_unique_acts))


bio_flows_correspondence = get_correspondence_bio_flows()
exc_codes = {}


@lru_cache
def fetch_exchange_code(name, ref, loc, unit):
    if (name, ref, loc, unit) not in exc_codes:
        code = str(uuid.uuid4().hex)
        exc_codes[(name, ref, loc, unit)] = code
    else:
        code = exc_codes[(name, ref, loc, unit)]

    return code


def get_act_dict_structure(ind, acts_ind, db_name) -> dict:
    name, ref, _, loc, unit, _ = acts_ind[ind]
    code = fetch_exchange_code(name, ref, loc, unit)

    return {
        "name": name,
        "reference product": ref,
        "unit": unit,
        "location": loc,
        "database": db_name,
        "code": code,
        "parameters": [],
        "exchanges": [],
    }


def correct_biosphere_flow(name, cat, unit, version):
    """
    Correct the biosphere flow name if it is outdated.
    """

    bio_dict = biosphere_flows_dictionary(version)

    try:

        if len(cat) > 1:
            main_cat = cat[0]
            sub_cat = cat[1]
        else:
            main_cat = cat[0]
            sub_cat = "unspecified"

        if (name, main_cat, sub_cat, unit) not in bio_dict:
            if bio_flows_correspondence.get(main_cat, {}).get(name, {}):
                name = bio_flows_correspondence[main_cat][name]
                return bio_dict[(name, main_cat, sub_cat, unit)]
    except:
        print(f"{name, cat, unit, version} not found in biosphere dictionary.")
        return None

    if (name, main_cat, sub_cat, unit) not in bio_dict:
        print(f"{name, cat, unit, version} not found in biosphere dictionary.")
        return None

    return bio_dict[(name, main_cat, sub_cat, unit)]


def get_exchange(ind, acts_ind, db_name, version, amount=1.0):
    name, ref, cat, loc, unit, flow_type = acts_ind[ind]
    _ = lambda x: x if x != 0 else 1.0
    return {
        "name": name,
        "product": ref,
        "unit": unit,
        "location": loc,
        "categories": cat,
        "type": flow_type,
        "amount": amount if flow_type != "production" else _(amount),
        "input": (
            ("biosphere3", correct_biosphere_flow(name, cat, unit, version))
            if flow_type == "biosphere"
            else (db_name, fetch_exchange_code(name, ref, loc, unit))
        ),
    }


def write_formatted_data(name, data, filepath):
    """
    Adapted from bw2io.export.csv
    :param name: name of the database
    :param data: data to write
    :param filepath: path to the file
    """

    sections = [
        "project parameters",
        "database",
        "database parameters",
        "activities",
        "activity parameters",
        "exchanges",
    ]

    result = []

    if "database" in sections:
        result.append(["Database", name])
        result.append([])

    if "activities" not in sections:
        return result

    for act in data:
        result.append(["Activity", act["name"]])
        result.append(["reference product", act["reference product"]])
        result.append(["unit", act["unit"]])
        result.append(["location", act["location"]])
        result.append(["comment", act.get("comment", "")])
        result.append(["source", act.get("source", "")])
        result.append(["parameters", act.get("parameters", [])])
        result.append([""])

        if "exchanges" in sections:
            result.append(["Exchanges"])
            if act.get("exchanges"):
                result.append(
                    [
                        "name",
                        "amount",
                        "unit",
                        "location",
                        "categories",
                        "type",
                        "product",
                    ]
                )
                for exc in act["exchanges"]:
                    result.append(
                        [
                            exc["name"],
                            exc["amount"],
                            exc["unit"],
                            exc.get("location"),
                            (
                                "::".join(list(exc.get("categories", [])))
                                if exc["type"] == "biosphere"
                                else None
                            ),
                            exc["type"],
                            exc.get("product"),
                        ]
                    )
        result.append([])

    with open(filepath, "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        for line in result:
            writer.writerow(line)

    return filepath


def build_datapackage(df, inventories, list_scenarios, ei_version, name):
    """
    Create and export a scenario datapackage.
    """

    # check that directory exists, otherwise create it
    Path(DIR_DATAPACKAGE_TEMP).mkdir(parents=True, exist_ok=True)
    df.to_csv(
        DIR_DATAPACKAGE_TEMP / "scenario_data.csv", index=False, encoding="utf-8-sig"
    )
    write_formatted_data(
        name=name, data=inventories, filepath=DIR_DATAPACKAGE_TEMP / "inventories.csv"
    )
    package = Package(base_path=str(DIR_DATAPACKAGE_TEMP))
    package.infer("**/*.csv")
    package.descriptor["name"] = name
    package.descriptor["title"] = name.capitalize()
    package.descriptor["description"] = (
        f"Data package generated by premise {__version__}."
    )
    package.descriptor["premise version"] = str(__version__)
    package.descriptor["dependencies"] = [
        {
            "name": "ecoinvent",
            "system model": "cut-off",
            "version": ei_version,
            "type": "source",
        },
        {
            "name": "biosphere3",
        },
    ]

    if len(list_scenarios[0]) == 3:
        package.descriptor["scenarios"] = [
            {
                "name": s,
                "description": f"Prospective db, "
                f"based on {s.split(' - ')[0].upper()}, "
                f"pathway {s.split(' - ')[1].upper()}, "
                f"for the year {s.split(' - ')[-1]}.",
            }
            for s in list_scenarios[1:]
        ]
    else:
        package.descriptor["scenarios"] = [
            {
                "name": s,
                "description": f"Prospective db, "
                f"based on {s.split(' - ')[0].upper()}, "
                f"pathway {s.split(' - ')[1].upper()}, "
                f"for the year {s.split(' - ')[2]}, and "
                f"external scenario {' '.join(s.split(' - ')[3:])}.",
            }
            for s in list_scenarios[1:]
        ]

    package.descriptor["keywords"] = [
        "ecoinvent",
        "scenario",
        "data package",
        "premise",
    ]
    package.descriptor["licenses"] = [
        {
            "id": "CC0-1.0",
            "title": "CC0 1.0",
            "url": "https://creativecommons.org/publicdomain/zero/1.0/",
        }
    ]
    package.commit()

    # check that directory exists, otherwise create it
    Path(DIR_DATAPACKAGE).mkdir(parents=True, exist_ok=True)

    # save the datapackage
    package.save(DIR_DATAPACKAGE / f"{name}.zip")

    print(f"Data package saved at {DIR_DATAPACKAGE / f'{name}.zip'}")


def generate_scenario_factor_file(
    origin_db: list,
    scenarios: list,
    db_name: str,
    biosphere_name: str,
    version: str,
    scenario_list: list = None,
):
    """
    Generate a scenario factor file from a list of databases
    :param origin_db: the original database
    :param scenarios: a list of databases
    :param db_name: the name of the database
    :param biosphere_name: the name of the biosphere database
    :param version: the version of ecoinvent
    :param scenario_list: a list of external scenarios
    """

    print("Building scenario factor file...")

    # create the dataframe
    df, new_db, list_unique_acts = generate_scenario_difference_file(
        origin_db=origin_db,
        scenarios=scenarios,
        db_name=db_name,
        biosphere_name=biosphere_name,
        version=version,
        scenario_list=scenario_list,
    )

    original = df["original"]
    original = original.replace(0, 1)
    df.loc[:, "original":] = df.loc[:, "original":].div(original, axis=0)

    # remove the column `original`
    df = df.drop(columns=["original"])

    # fetch a list of activities not present in original_db
    list_original_acts = get_list_unique_acts([{"database": origin_db}])

    new_acts_list = list(set(list_unique_acts) - set(list_original_acts))

    print(f"Number of new activities: {len(new_acts_list)}")

    # turn new_acts_list into a dictionary
    new_acts_dict = {v: k for k, v in dict(enumerate(new_acts_list)).items()}

    # fetch the additional activities from new_db
    extra_acts = [
        dataset
        for dataset in new_db
        if (
            dataset["name"],
            dataset.get("reference product"),
            None,
            dataset.get("location"),
            dataset["unit"],
            "production",
        )
        in new_acts_dict
    ]

    return df, extra_acts


def generate_new_activities(args):
    k, v, acts_ind, db_name, version, dict_meta, m = args
    act = get_act_dict_structure(k, acts_ind, db_name)
    meta_id = tuple(list(acts_ind[k])[:-1])
    act.update(dict_meta[meta_id])

    act["exchanges"].extend(
        get_exchange(i, acts_ind, db_name, version, amount=m[i, k, 0]) for i in v
    )
    return act


def generate_scenario_difference_file(
    db_name, origin_db, scenarios, version, scenario_list, biosphere_name
) -> tuple[DataFrame, list[dict], set[Any]]:
    """
    Generate a scenario difference file for a given list of databases
    :param db_name: name of the new database
    :param origin_db: the original database
    :param scenarios: list of databases
    """

    bio_dict = biosphere_flows_dictionary(version)

    exc_codes.update(
        {
            (a["name"], a["reference product"], a["location"], a["unit"]): a["code"]
            for a in origin_db
        }
    )
    # Turn list into set for O(1) membership tests
    list_acts = set(get_list_unique_acts([{"database": origin_db}] + scenarios))

    acts_ind = dict(enumerate(list_acts))
    acts_ind_rev = {v: k for k, v in acts_ind.items()}

    list_scenarios = ["original"] + scenario_list

    list_dbs = [origin_db] + [a["database"] for a in scenarios]

    matrices = {
        a: nsp.lil_matrix((len(list_acts), len(list_acts)))
        for a, _ in enumerate(list_scenarios)
    }

    # Use defaultdict to avoid key errors
    dict_meta = defaultdict(dict)

    for db in list_dbs:
        for a in db:
            key = (a["name"], a["reference product"], None, a["location"], a["unit"])
            dict_meta[key] = {
                b: c
                for b, c in a.items()
                if b
                not in [
                    "exchanges",
                    "code",
                    "name",
                    "reference product",
                    "location",
                    "unit",
                    "database",
                ]
            }

    for i, db in enumerate(list_dbs):
        for ds in db:
            c = (
                ds["name"],
                ds.get("reference product"),
                ds.get("categories"),
                ds.get("location"),
                ds["unit"],
                "production",
            )
            for exc in ds["exchanges"]:
                s = (
                    exc["name"],
                    exc.get("product"),
                    exc.get("categories"),
                    exc.get("location"),
                    exc["unit"],
                    exc["type"],
                )
                matrices[i][acts_ind_rev[s], acts_ind_rev[c]] += exc["amount"]

    m = sparse.stack([sparse.COO(x) for x in matrices.values()], axis=-1)
    inds = sparse.argwhere(m.sum(-1).T != 0)
    inds = list(map(tuple, inds))

    dataframe_rows = []

    inds_d = defaultdict(list)
    for ind in inds:
        inds_d[ind[0]].append(ind[1])

    with Pool(processes=mp.cpu_count()) as pool:
        new_db = pool.map(
            generate_new_activities,
            [
                (k, v, acts_ind, db_name, version, dict_meta, m)
                for k, v in inds_d.items()
            ],
        )

    inds_std = sparse.argwhere((m[..., 1:] == m[..., 0, None]).all(axis=-1).T == False)

    for i in inds_std:
        c_name, c_ref, c_cat, c_loc, c_unit, _ = acts_ind[i[0]]
        s_name, s_ref, s_cat, s_loc, s_unit, s_type = acts_ind[i[1]]

        database_name = db_name
        exc_key_supplier = None

        if s_type == "biosphere":
            database_name = "biosphere3"

            key_exc = (
                s_name,
                s_cat[0],
                s_cat[1] if len(s_cat) > 1 else "unspecified",
                s_unit,
            )

            if key_exc in bio_dict:
                exc_key_supplier = (
                    database_name,
                    bio_dict[key_exc],
                )
            else:
                exc_key_supplier = (
                    database_name,
                    bio_dict[
                        bio_flows_correspondence.get(s_cat[0], {}).get(s_name, s_name),
                        s_cat[0],
                        s_cat[1] if len(s_cat) > 1 else "unspecified",
                        s_unit,
                    ],
                )

        row = [
            s_name,
            s_ref,
            s_loc,
            s_cat,
            database_name,
            exc_key_supplier,
            s_unit,
            c_name,
            c_ref,
            c_loc,
            c_cat,
            c_unit,
            db_name,
            None,
            s_type,
        ]

        row.extend(m[i[1], i[0], :])

        dataframe_rows.append(row)

    columns = [
        "from activity name",
        "from reference product",
        "from location",
        "from categories",
        "from database",
        "from key",
        "from unit",
        "to activity name",
        "to reference product",
        "to location",
        "to categories",
        "to unit",
        "to database",
        "to key",
        "flow type",
    ] + list_scenarios

    df = pd.DataFrame(dataframe_rows, columns=columns)

    df["to categories"] = None
    df = df.replace({"None": None, np.nan: None})
    df.loc[
        df["flow type"] == "biosphere", ["from reference product", "from location"]
    ] = None
    df.loc[df["flow type"].isin(["technosphere", "production"]), "from categories"] = (
        None
    )
    df.loc[df["flow type"] == "production", list_scenarios] = 1.0

    df.loc[df["flow type"] == "biosphere", "from database"] = biosphere_name

    # update the tuples in `from key` to make sure the first element
    # is the biosphere database name
    df.loc[df["flow type"] == "biosphere", "from key"] = df.loc[
        df["flow type"] == "biosphere", "from key"
    ].map(lambda x: (biosphere_name, x[1]))

    new_db, df = find_technosphere_keys(new_db, df)

    # return the dataframe and the new db
    return df, new_db, list_acts


def find_technosphere_keys(db, df):
    # erase keys for technosphere and production exchanges
    df.loc[df["flow type"].isin(["technosphere", "production"]), "from key"] = None
    df.loc[df["flow type"].isin(["technosphere", "production"]), "to key"] = None
    df.loc[df["flow type"] == "biosphere", "to key"] = None

    # reset all codes
    db = reset_all_codes(db)

    # create a dictionary of all activities
    dict_act = {
        (a["name"], a["reference product"], a["location"]): (a["database"], a["code"])
        for a in db
    }

    # iterate through df
    # and fill "from key" and "to key" columns
    # if None

    df.loc[df["from key"].isnull(), "from key"] = pd.Series(
        list(
            zip(
                df["from activity name"],
                df["from reference product"],
                df["from location"],
            )
        )
    ).map(dict_act)

    df.loc[df["to key"].isnull(), "to key"] = pd.Series(
        list(zip(df["to activity name"], df["to reference product"], df["to location"]))
    ).map(dict_act)

    return db, df


def generate_superstructure_db(
    origin_db,
    scenarios,
    db_name,
    biosphere_name,
    filepath,
    version,
    scenario_list,
    file_format="excel",
) -> List[dict]:
    """
    Build a superstructure database from a list of databases
    :param origin_db: the original database
    :param scenarios: a list of modified databases
    :param db_name: the name of the new database
    :param filepath: the filepath of the new database
    :param version: the version of the new database
    :param file_format: the format of the scenario difference file. Can be "excel", "csv" or "feather".
    :return: a superstructure database
    """

    print("Building superstructure database...")

    # create the dataframe
    df, new_db, _ = generate_scenario_difference_file(
        origin_db=origin_db,
        scenarios=scenarios,
        db_name=db_name,
        biosphere_name=biosphere_name,
        version=version,
        scenario_list=scenario_list,
    )

    # remove unneeded columns "to unit"
    df = df.drop(columns=["to unit"])

    # rename column "from unit" to "unit"
    df = df.rename(columns={"from unit": "unit"})

    # remove the column `original`
    df = df.drop(columns=["original"])
    if "unit" in df.columns:
        df = df.drop(columns=["unit"])

    if filepath is not None:
        filepath = Path(filepath)
    else:
        filepath = Path.cwd() / "export" / "scenario diff files"

    if not os.path.exists(filepath):
        os.makedirs(filepath)

    # Drop duplicate rows
    # should not be any, but just in case
    before = len(df)
    df = df.drop_duplicates()
    after = len(df)
    print(f"Dropped {before - after} duplicate(s).")

    # if df is longer than the row limit of Excel,
    # the export to Excel is not an option
    if len(df) > 1048576:
        file_format = "csv"
        print(
            "The scenario difference file is too long to be exported to Excel. Exporting to CSV instead."
        )

    if file_format == "excel":
        filepath_sdf = filepath / f"scenario_diff_{db_name}.xlsx"
        df.to_excel(filepath_sdf, index=False)
    elif file_format == "csv":
        filepath_sdf = filepath / f"scenario_diff_{db_name}.csv"
        df.to_csv(filepath_sdf, index=False, sep=";", encoding="utf-8-sig")
    elif file_format == "feather":
        filepath_sdf = filepath / f"scenario_diff_{db_name}.feather"
        df.to_feather(filepath_sdf)
    else:
        raise ValueError(f"Unknown format {file_format}")

    print(f"Scenario difference file exported to {filepath}!")

    return new_db


def check_geographical_linking(scenario, original_database):

    # geo = Geomap(scenario["model"])

    index = scenario.get("index") or {}
    database = scenario["database"]
    original_datasets = [
        (a["name"], a["reference product"], a["location"]) for a in original_database
    ]

    datasets_to_check = [
        ds
        for ds in database
        if (ds["name"], ds["reference product"], ds["location"])
        not in original_datasets
    ]

    FORBIDDEN = ["import", "mix", "transport"]

    for ds in datasets_to_check:
        if ds["location"] not in ["GLO", "RoW", "World"]:
            if "market" not in ds["name"]:
                for exc in ds["exchanges"]:
                    if exc["type"] == "technosphere":
                        if not any(x in exc["name"] for x in FORBIDDEN) and not any(
                            x in ds["name"] for x in FORBIDDEN
                        ):
                            if exc["location"] != ds["location"]:
                                # check if exchange from the same location as the dataset is available
                                key = (exc["name"], exc["product"])
                                if ds["location"] in [
                                    k["location"] for k in index.get(key, [])
                                ]:
                                    # if ds["location"] not in geo.iam_to_ecoinvent_location(exc["location"]):
                                    if (exc["name"], exc["product"]) != (
                                        ds["name"],
                                        ds["reference product"],
                                    ):
                                        # there is a better match available
                                        exc["location"] = ds["location"]

    return scenario


def prepare_db_for_export(
    scenario, name, original_database, keep_uncertainty_data=False, biosphere_name=None
):
    """
    Prepare a database for export.
    """

    # ensuring that all geographically appropriate exchanges are present
    scenario = check_geographical_linking(scenario, original_database)

    # validate the database
    validator = BaseDatasetValidator(
        model=scenario["model"],
        scenario=scenario["pathway"],
        year=scenario["year"],
        regions=scenario["iam data"].regions,
        original_database=original_database,
        database=scenario["database"],
        db_name=name,
        keep_uncertainty_data=keep_uncertainty_data,
        biosphere_name=biosphere_name,
    )
    validator.run_all_checks()

    return validator.database


def _prepare_database(
    scenario, db_name, original_database, keep_uncertainty_data, biosphere_name
):

    scenario["database"] = prepare_db_for_export(
        scenario,
        name=db_name,
        original_database=original_database,
        keep_uncertainty_data=keep_uncertainty_data,
        biosphere_name=biosphere_name,
    )

    return scenario


class Export:
    """
    Class that exports the transformed data into matrices:
    * A matrix: contains products exchanges
    * B matrix: contains exchanges activities and the biosphere

    The A and B matrices are exported as csv files in a sparse representation (only non-zero values are listed), like so:
    - index row, index column, value of exchange

    Dictionaries to map row numbers to activities and products names are also exported.

    :ivar db: transformed database
    :vartype database: dict
    :ivar scenario: name of an IAM pathway
    :vartype pathway: str
    :ivar year: year of scenario
    :vartype year: int

    """

    def __init__(
        self,
        scenario: dict = None,
        filepath: Path = None,
        version: str = None,
    ):
        self.db = scenario["database"]
        self.model = scenario["model"]
        self.scenario = scenario["pathway"]
        self.year = scenario["year"]
        self.filepath = filepath
        self.version = version
        self.bio_codes = self.rev_index(
            create_codes_index_of_biosphere_flows_matrix(self.version)
        )
        self.bio_dict = biosphere_flows_dictionary(self.version)
        self.unmatched_category_flows = []

    def create_A_matrix_coordinates(self) -> list:
        """
        Create the coordinates of the A matrix.

        :param export_uncertainty: if True, export uncertainty data
        """
        index_A = create_index_of_A_matrix(self.db)
        list_exchanges = []

        try:
            for ds in self.db:
                for exc in ds["exchanges"]:

                    if exc["type"] == "production":
                        row = [
                            index_A[
                                (
                                    ds["name"],
                                    ds["reference product"],
                                    ds["unit"],
                                    ds["location"],
                                )
                            ],
                            index_A[
                                (
                                    exc["name"],
                                    exc["product"],
                                    exc["unit"],
                                    exc["location"],
                                )
                            ],
                            exc["amount"],
                            0,
                            exc["amount"],
                            None,
                            None,
                            None,
                            None,
                            0,
                            0,
                        ]
                        list_exchanges.append(row)

                    if exc["type"] == "technosphere":
                        row = [
                            index_A[
                                (
                                    ds["name"],
                                    ds["reference product"],
                                    ds["unit"],
                                    ds["location"],
                                )
                            ],
                            index_A[
                                (
                                    exc["name"],
                                    exc["product"],
                                    exc["unit"],
                                    exc["location"],
                                )
                            ],
                            exc["amount"],
                            exc.get("uncertainty type", 0),
                            (
                                exc["loc"]
                                if exc.get("uncertainty type", 0) not in [0, 1]
                                else exc["amount"]
                            ),
                            exc.get("scale"),
                            exc.get("shape"),
                            exc.get("minimum"),
                            exc.get("maximum"),
                            1 if exc["amount"] < 0 else 0,
                            1,
                        ]
                        list_exchanges.append(row)

        except KeyError:
            print(f"KeyError for {exc} in {ds['name']}")

        return list_exchanges

    def create_B_matrix_coordinates(self):
        index_B = create_index_of_biosphere_flows_matrix(self.version)
        rev_index_B = self.create_rev_index_of_B_matrix(self.version)
        index_A = create_index_of_A_matrix(self.db)
        list_rows = []

        for ds in self.db:
            for exc in ds["exchanges"]:
                if exc["type"] == "biosphere":
                    try:
                        code = exc["input"][1]
                        lookup = rev_index_B[code]
                        ind_B = index_B[lookup]

                        row = [
                            index_A[
                                (
                                    ds["name"],
                                    ds["reference product"],
                                    ds["unit"],
                                    ds["location"],
                                )
                            ],
                            ind_B,
                            exc["amount"],
                            exc.get("uncertainty type", 0),
                            (
                                exc["loc"]
                                if exc.get("uncertainty type", 0) not in [0, 1]
                                else exc["amount"]
                            ),
                            exc.get("scale"),
                            exc.get("shape"),
                            exc.get("minimum"),
                            exc.get("maximum"),
                            1 if exc["amount"] < 0 else 0,
                            1,
                        ]
                    except KeyError:
                        print(
                            "Cannot find the biosphere flow",
                            exc["name"],
                            exc["categories"],
                        )
                        row = ()
                    list_rows.append(row)
        return list_rows

    def export_db_to_matrices(self):
        """
        Export the database to A and B matrices.

        :param export_uncertainty: if True, export uncertainty data
        """

        if not os.path.exists(self.filepath):
            os.makedirs(self.filepath)

        # Export A matrix
        rows = self.create_A_matrix_coordinates()

        with open(self.filepath / "A_matrix.csv", "w", encoding="utf-8") as file:
            writer = csv.writer(
                file,
                delimiter=";",
                lineterminator="\n",
            )
            writer.writerow(
                [
                    "index of activity",
                    "index of product",
                    "value",
                    "uncertainty type",
                    "loc",
                    "scale",
                    "shape",
                    "minimum",
                    "maximum",
                    "negative",
                    "flip",
                ]
            )
            for row in rows:
                writer.writerow(row)

        # Export A index
        with open(self.filepath / "A_matrix_index.csv", "w", encoding="utf-8") as file:
            writer = csv.writer(
                file,
                delimiter=";",
                lineterminator="\n",
            )
            index_A = create_index_of_A_matrix(self.db)
            for d in index_A:
                data = list(d) + [index_A[d]]
                writer.writerow(data)

        index_B = create_index_of_biosphere_flows_matrix(self.version)

        # Export B matrix
        with open(self.filepath / "B_matrix.csv", "w", encoding="utf-8") as file:
            writer = csv.writer(
                file,
                delimiter=";",
                lineterminator="\n",
            )
            writer.writerow(
                [
                    "index of activity",
                    "index of biosphere flow",
                    "value",
                    "uncertainty type",
                    "loc",
                    "scale",
                    "shape",
                    "minimum",
                    "maximum",
                    "negative",
                    "flip",
                ]
            )
            rows = self.create_B_matrix_coordinates()
            for row in rows:
                writer.writerow(row)

        # Export B index
        with open(self.filepath / "B_matrix_index.csv", "w", encoding="utf-8") as file:
            writer = csv.writer(
                file,
                delimiter=";",
                lineterminator="\n",
            )
            for d in index_B:
                data = list(d) + [index_B[d]]
                writer.writerow(data)

        print(f"Matrices saved in {self.filepath}.")

    @staticmethod
    def create_rev_index_of_B_matrix(version):
        return {v: k for k, v in biosphere_flows_dictionary(version).items()}

    def get_category_of_exchange(self):
        """
        This function returns a dictionnary with (name, reference product) as keys,
        and {'main category': 'xxxxxx', 'category': 'yyyyyy'} as values.
        This is useful for Simapro export, to categorize datasets into a tree structure.
        :return: dict
        """

        dict_classifications = load_simapro_categories()

        dict_categories = {}

        for ds in self.db:
            if (ds["name"], ds["reference product"]) not in dict_categories:
                main_category, category = (None, None)

                if "classifications" in ds:
                    if len(ds["classifications"]) > 0:
                        for x in ds["classifications"]:
                            if x[0] == "ISIC rev.4 ecoinvent":
                                try:
                                    key = [
                                        int(s)
                                        for s in re.findall(r"\d+", x[1].split(":")[0])
                                    ][0]
                                    main_category = dict_classifications[key][
                                        "category 1"
                                    ]

                                except KeyError:
                                    try:
                                        key = x[1].split(":")[0].strip()
                                        main_category = dict_classifications[key][
                                            "category 1"
                                        ]

                                    except KeyError:
                                        continue

                                if dict_classifications[key]["category 3"] != "":
                                    category = (
                                        dict_classifications[key]["category 2"]
                                        + "\ ".strip()
                                        + dict_classifications[key]["category 3"]
                                    )
                                else:
                                    category = dict_classifications[key]["category 2"]

                        if not main_category:
                            for x in ds["classifications"]:
                                if x[0] == "CPC":
                                    main_category = dict_classifications[
                                        x[1].split(":")[0].strip()
                                    ]["category 1"]

                                    if (
                                        dict_classifications[
                                            x[1].split(":")[0].strip()
                                        ]["category 3"]
                                        != ""
                                    ):
                                        category = (
                                            dict_classifications[
                                                x[1].split(":")[0].strip()
                                            ]["category 2"]
                                            + "\ ".strip()
                                            + dict_classifications[
                                                x[1].split(":")[0].strip()
                                            ]["category 3"]
                                        )
                                    else:
                                        category = dict_classifications[
                                            x[1].split(":")[0].strip()
                                        ]["category 2"]

                if not main_category or main_category == "":
                    main_category = "material"
                    category = "Others"

                dict_categories[(ds["name"], ds["reference product"])] = {
                    "main category": main_category,
                    "category": category,
                }

        return dict_categories

    def export_db_to_simapro(self, olca_compartments=False):
        if not os.path.exists(self.filepath):
            os.makedirs(self.filepath)

        dict_bio = get_simapro_biosphere_dictionnary()

        headers = [
            "{SimaPro 9.1.1.7}",
            "{processes}",
            "{Project: premise import" + f"{datetime.today():%d.%m.%Y}" + "}",
            "{CSV Format version: 9.0.0}",
            "{CSV separator: Semicolon}",
            "{Decimal separator: .}",
            "{Date separator: .}",
            "{Short date format: dd.MM.yyyy}",
            "{Export platform IDs: No}",
            "{Skip empty fields: No}",
            "{Convert expressions to constants: No}",
            "{Related objects(system descriptions, substances, units, etc.): Yes}",
            "{Include sub product stages and processes: Yes}",
        ]

        fields = [
            "Process",
            "Category type",
            "Type",
            "Process name",
            "Time Period",
            "Geography",
            "Technology",
            "Representativeness",
            "Waste treatment allocation",
            "Cut off rules",
            "Capital goods",
            "Date",
            "Boundary with nature",
            "Infrastructure",
            "Record",
            "Generator",
            "Literature references",
            "External documents",
            "Comment",
            "Collection method",
            "Data treatment",
            "Verification",
            "System description",
            "Allocation rules",
            "Products",
            "Waste treatment",
            "Materials/fuels",
            "Resources",
            "Emissions to air",
            "Emissions to water",
            "Emissions to soil",
            "Non material emission",
            "Social issues",
            "Economic issues",
            "Waste to treatment",
            "End",
        ]

        # mapping between BW2 and Simapro units
        simapro_units = get_simapro_units()

        # mapping between BW2 and Simapro sub-compartments
        if olca_compartments:
            simapro_subs = {}
        else:
            simapro_subs = get_simapro_compartments()

        filename = f"simapro_export_{self.model}_{self.scenario}_{self.year}.csv"

        dict_cat_simapro = get_simapro_category_of_exchange()
        dict_cat = self.get_category_of_exchange()
        dict_refs = load_references()

        unlinked_biosphere_flows = []

        with open(
            Path(self.filepath) / filename, "w", newline="", encoding="latin1"
        ) as csvFile:
            writer = csv.writer(csvFile, delimiter=";")
            for item in headers:
                writer.writerow([item])
            writer.writerow([])

            for ds in self.db:
                key = (ds["name"].lower(), ds["reference product"].lower())

                try:
                    main_category, sub_category = (
                        dict_cat_simapro[key]["category"],
                        dict_cat_simapro[key]["sub_category"],
                    )
                except KeyError:
                    main_category, sub_category = ("material", "Others\Transformation")
                    self.unmatched_category_flows.append(key)

                for item in fields:
                    if (
                        main_category.lower() == "waste treatment"
                        and item == "Products"
                    ):
                        continue

                    if main_category.lower() != "waste treatment" and item in (
                        "Waste treatment",
                        "Waste treatment allocation",
                    ):
                        continue

                    writer.writerow([item])

                    if item == "Process name":
                        name = f"{ds['reference product']} {{{ds.get('location', 'GLO')}}}| {ds['name']} | Cut-off, U"

                        writer.writerow([name])

                    if item == "Type":
                        writer.writerow(["Unit process"])

                    if item == "Category type":
                        writer.writerow([main_category])

                    if item == "Generator":
                        writer.writerow(["premise " + str(__version__)])

                    if item == "Geography":
                        writer.writerow([ds["location"]])

                    if item == "Date":
                        writer.writerow([f"{datetime.today():%d.%m.%Y}"])

                    if item == "Comment":
                        if ds["name"] in dict_refs:
                            string = re.sub(
                                "[^a-zA-Z0-9 .,]", "", dict_refs[ds["name"]]["source"]
                            )

                            if dict_refs[ds["name"]]["description"] != "":
                                string += " " + re.sub(
                                    "[^a-zA-Z0-9 .,]",
                                    "",
                                    dict_refs[ds["name"]]["description"],
                                )

                            writer.writerow([string])
                        else:
                            if "comment" in ds:
                                string = re.sub("[^a-zA-Z0-9 .,]", "", ds["comment"])
                                writer.writerow([string])

                    if item in (
                        "Cut off rules",
                        "Capital goods",
                        "Technology",
                        "Representativeness",
                        "Waste treatment allocation",
                        "Boundary with nature",
                        "Allocation rules",
                        "Collection method",
                        "Verification",
                        "Time Period",
                        "Record",
                    ):
                        writer.writerow(["Unspecified"])
                    if item == "System description":
                        writer.writerow(["Ecoinvent v3"])
                    if item == "Infrastructure":
                        writer.writerow(["Yes"])
                    if item == "External documents":
                        writer.writerow(
                            [
                                "https://premise.readthedocs.io/en/latest/introduction.html"
                            ]
                        )
                    if item in ("Waste treatment", "Products"):
                        for e in ds["exchanges"]:
                            if e["type"] == "production":
                                name = f"{e['product']} {{{e.get('location', 'GLO')}}}| {e['name']} | Cut-off, U"

                                if item == "Waste treatment":
                                    writer.writerow(
                                        [
                                            name,
                                            simapro_units[e["unit"]],
                                            1.0,
                                            "not defined",
                                            sub_category,
                                        ]
                                    )

                                else:
                                    writer.writerow(
                                        [
                                            name,
                                            simapro_units[e["unit"]],
                                            1.0,
                                            "100%",
                                            "not defined",
                                            sub_category,
                                        ]
                                    )
                                e["used"] = True
                    if item == "Materials/fuels":
                        for e in ds["exchanges"]:
                            if e["type"] == "technosphere":
                                key = (e["name"].lower(), e["product"].lower())
                                if key in dict_cat_simapro:
                                    exc_cat = dict_cat_simapro.get(key)["category"]
                                else:
                                    exc_cat = "material"

                                if exc_cat != "waste treatment":
                                    name = f"{e['product']} {{{e.get('location', 'GLO')}}}| {e['name']} | Cut-off, U"

                                    writer.writerow(
                                        [
                                            name,
                                            simapro_units[e["unit"]],
                                            f"{e['amount']:.3E}",
                                            "undefined",
                                            0,
                                            0,
                                            0,
                                        ]
                                    )
                                    e["used"] = True
                    if item == "Resources":
                        for e in ds["exchanges"]:
                            if (
                                e["type"] == "biosphere"
                                and e["categories"][0] == "natural resource"
                            ):
                                if e["name"] not in dict_bio:
                                    unlinked_biosphere_flows.append(
                                        [e["name"], " - ", e["categories"][0]]
                                    )

                                if len(e["categories"]) > 1:
                                    if e["categories"][1] != "fossil well":
                                        sub_compartment = simapro_subs.get(
                                            e["categories"][1], e["categories"][1]
                                        )
                                    else:
                                        sub_compartment = ""

                                writer.writerow(
                                    [
                                        dict_bio.get(e["name"], e["name"]),
                                        sub_compartment,
                                        simapro_units[e["unit"]],
                                        f"{e['amount']:.3E}",
                                        "undefined",
                                        0,
                                        0,
                                        0,
                                    ]
                                )
                                e["used"] = True
                    if item == "Emissions to air":
                        for e in ds["exchanges"]:
                            if e["type"] == "biosphere" and e["categories"][0] == "air":
                                if len(e["categories"]) > 1:
                                    sub_compartment = simapro_subs.get(
                                        e["categories"][1], e["categories"][1]
                                    )
                                else:
                                    sub_compartment = ""

                                if e["name"].lower() == "water":
                                    e["unit"] = "kilogram"
                                    # going from cubic meters to kilograms
                                    e["amount"] *= 1000

                                if e["name"] not in dict_bio:
                                    unlinked_biosphere_flows.append(
                                        [e["name"], " - ", e["categories"][0]]
                                    )

                                writer.writerow(
                                    [
                                        dict_bio.get(e["name"], e["name"]),
                                        sub_compartment,
                                        simapro_units[e["unit"]],
                                        f"{e['amount']:.3E}",
                                        "undefined",
                                        0,
                                        0,
                                        0,
                                    ]
                                )
                                e["used"] = True
                    if item == "Emissions to water":
                        for e in ds["exchanges"]:
                            if (
                                e["type"] == "biosphere"
                                and e["categories"][0] == "water"
                            ):
                                if len(e["categories"]) > 1:
                                    sub_compartment = simapro_subs.get(
                                        e["categories"][1], e["categories"][1]
                                    )
                                else:
                                    sub_compartment = ""

                                if e["name"].lower() == "water":
                                    e["unit"] = "kilogram"
                                    e["amount"] /= 1000

                                if e["name"] not in dict_bio:
                                    unlinked_biosphere_flows.append(
                                        [e["name"], " - ", e["categories"][0]]
                                    )

                                writer.writerow(
                                    [
                                        dict_bio.get(e["name"], e["name"]),
                                        sub_compartment,
                                        simapro_units[e["unit"]],
                                        f"{e['amount']:.3E}",
                                        "undefined",
                                        0,
                                        0,
                                        0,
                                    ]
                                )
                                e["used"] = True
                    if item == "Emissions to soil":
                        for e in ds["exchanges"]:
                            if (
                                e["type"] == "biosphere"
                                and e["categories"][0] == "soil"
                            ):
                                if len(e["categories"]) > 1:
                                    sub_compartment = simapro_subs.get(
                                        e["categories"][1], e["categories"][1]
                                    )
                                else:
                                    sub_compartment = ""

                                if e["name"] not in dict_bio:
                                    unlinked_biosphere_flows.append(
                                        [e["name"], " - ", e["categories"][0]]
                                    )

                                writer.writerow(
                                    [
                                        dict_bio.get(e["name"], e["name"]),
                                        sub_compartment,
                                        simapro_units[e["unit"]],
                                        f"{e['amount']:.3E}",
                                        "undefined",
                                        0,
                                        0,
                                        0,
                                    ]
                                )
                                e["used"] = True
                    if item == "Waste to treatment":
                        for e in ds["exchanges"]:
                            if e["type"] == "technosphere":
                                key = (e["name"].lower(), e["product"].lower())
                                if key in dict_cat_simapro:
                                    exc_cat = dict_cat_simapro.get(key)["category"]
                                else:
                                    exc_cat = "material"

                                if exc_cat == "waste treatment":
                                    name = f"{e['product']} {{{e.get('location', 'GLO')}}}| {e['name']} | Cut-off, U"

                                    writer.writerow(
                                        [
                                            name,
                                            simapro_units[e["unit"]],
                                            f"{e['amount'] * -1:.3E}",
                                            "undefined",
                                            0,
                                            0,
                                            0,
                                        ]
                                    )
                                    e["used"] = True

                    writer.writerow([])

            # System description
            writer.writerow(["System description"])
            writer.writerow([])
            writer.writerow(["Name"])
            writer.writerow(["Ecoinvent v3"])
            writer.writerow([])
            writer.writerow(["Category"])
            writer.writerow(["Others"])
            writer.writerow([])
            writer.writerow(["Description"])
            writer.writerow([""])
            writer.writerow([])
            writer.writerow(["Cut-off rules"])
            writer.writerow([""])
            writer.writerow([])
            writer.writerow(["Energy model"])
            writer.writerow([])
            writer.writerow([])
            writer.writerow(["Transport model"])
            writer.writerow([])
            writer.writerow([])
            writer.writerow(["Allocation rules"])
            writer.writerow([])
            writer.writerow(["End"])
            writer.writerow([])

        csvFile.close()

        # check that all exchanges have been used
        unused_exchanges = []
        for ds in self.db:
            for e in ds["exchanges"]:
                if "used" not in e:
                    if [
                        e["name"][:40],
                        e.get("product", "")[:40],
                        e.get("categories"),
                        e.get("location"),
                    ] not in unused_exchanges:
                        unused_exchanges.append(
                            [
                                e["name"][:40],
                                e.get("product", "")[:40],
                                e.get("categories"),
                                e.get("location"),
                            ]
                        )

        if len(unused_exchanges) > 0:
            print("The following exchanges have not been used in the Simapro export:")
            # make prettytable
            x = PrettyTable()
            x.field_names = ["Name", "Product", "Categories", "Location"]
            for i in unused_exchanges:
                x.add_row(i)
            print(x)

        if len(self.unmatched_category_flows) > 0:
            print(f"{len(self.unmatched_category_flows)} unmatched flow categories.")

        print(f"Simapro CSV file saved in {self.filepath}.")

    def rev_index(self, inds):
        return {v: k for k, v in inds.items()}

    def get_bio_code(self, idx):
        return self.bio_codes[idx]
