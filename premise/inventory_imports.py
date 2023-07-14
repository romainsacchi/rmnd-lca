"""
Contains class and methods to imports inventories from ecoinvent, premise,
and those provided by the user.
"""

import csv
import itertools
import sys
import uuid
from functools import lru_cache
from pathlib import Path
from typing import Dict, List, Union

import bw2io
import numpy as np
import requests
import yaml
from bw2io import CSVImporter, ExcelImporter, Migration
from prettytable import PrettyTable
from wurst import searching as ws

from . import DATA_DIR, INVENTORY_DIR
from .clean_datasets import remove_categories, remove_uncertainty
from .data_collection import get_delimiter
from .geomap import Geomap

FILEPATH_MIGRATION_MAP = INVENTORY_DIR / "migration_map.csv"
FILEPATH_CONSEQUENTIAL_BLACKLIST = DATA_DIR / "consequential" / "blacklist.yaml"
CORRESPONDENCE_BIO_FLOWS = (
    DATA_DIR / "utils" / "export" / "correspondence_biosphere_flows.yaml"
)


def get_correspondence_bio_flows():
    """
    Mapping between ei39 and ei<39 biosphere flows.
    """

    with open(CORRESPONDENCE_BIO_FLOWS, "r", encoding="utf-8") as stream:
        flows = yaml.safe_load(stream)

    return flows


@lru_cache
def get_biosphere_code(version) -> dict:
    """
    Retrieve a dictionary with biosphere flow names and uuid codes.
    :returns: dictionary with biosphere flow names as keys and uuid codes as values

    """
    if version == "3.9":
        fp = DATA_DIR / "utils" / "export" / "flows_biosphere_39.csv"
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
            csv_dict[(row[0], row[1], row[2], row[3])] = row[4]

    return csv_dict


def get_consequential_blacklist():
    with open(FILEPATH_CONSEQUENTIAL_BLACKLIST, "r", encoding="utf-8") as stream:
        flows = yaml.safe_load(stream)

    return flows


@lru_cache
def generate_migration_maps(origin: str, destination: str) -> Dict[str, list]:
    """
    Generate mapping for ecoinvent datasets across different database versions.
    :param origin: ecoinvent database version to find equivalence from (e.g., "3.6")
    :param destination: ecoinvent database version to find equivalence for (e.g., "3.8")
    :return: a migration map dicitonary for bw2io
    """

    response = {"fields": ["name", "reference product", "location"], "data": []}

    with open(FILEPATH_MIGRATION_MAP, "r", encoding="utf-8") as read_obj:
        csv_reader = csv.reader(
            read_obj,
            delimiter=get_delimiter(filepath=FILEPATH_MIGRATION_MAP),
        )
        next(csv_reader)
        for row in csv_reader:
            if row[0] == origin and row[1] == destination:
                data = {}
                if row[5] != "":
                    data["name"] = row[5]
                if row[6] != "":
                    data["reference product"] = row[6]
                if row[7] != "":
                    data["location"] = row[7]
                response["data"].append(((row[2], row[3], row[4]), data))

            if row[0] == destination and row[1] == origin:
                data = {}
                if row[2] != "":
                    data["name"] = row[2]
                if row[3] != "":
                    data["reference product"] = row[3]
                if row[4] != "":
                    data["location"] = row[4]
                response["data"].append(((row[5], row[6], row[7]), data))

    return response


def check_for_duplicate_datasets(data: List[dict]) -> List[dict]:
    """Check whether there are duplicate datasets in the inventory to import."""
    datasets = [(ds["name"], ds["reference product"], ds["location"]) for ds in data]
    duplicates = [
        item
        for item, count in itertools.groupby(sorted(datasets))
        if len(list(count)) > 1
    ]
    if duplicates:
        print("Duplicate datasets found (they need to be removed):")
        # print them using prettytable
        table = PrettyTable()
        table.field_names = ["Name", "Reference product", "Location"]
        for duplicate in duplicates:
            table.add_row(duplicate)
        print(table)

        # remove duplicates
        duplicates_added = []
        for ds in data:
            if (ds["name"], ds["reference product"], ds["location"]) in duplicates:
                if (
                    ds["name"],
                    ds["reference product"],
                    ds["location"],
                ) not in duplicates_added:
                    duplicates_added.append(
                        (ds["name"], ds["reference product"], ds["location"])
                    )
                else:
                    data.remove(ds)

    return data


def check_for_datasets_compliance_with_consequential_database(
    datasets: List[dict], blacklist: List[dict]
):
    """
    Check whether the datasets to import are compliant with the consequential database.

    :param datasets: list of datasets to import
    :param blacklist: list of datasets that are not in the consequential database
    :return: list of datasets that are compliant with the consequential database

    """
    # if system model is `consequential`` there is a
    # number of datasets we do not want to import

    tuples_of_blacklisted_datasets = [
        (i["name"], i["reference product"], i["unit"]) for i in blacklist
    ]

    datasets = [
        d
        for d in datasets
        if (d["name"], d["reference product"], d["unit"])
        not in tuples_of_blacklisted_datasets
    ]

    # also, we want to change exchanges that do not
    # exist in the consequential LCA database
    # and change them for the consequential equivalent

    for ds in datasets:
        for exchange in ds["exchanges"]:
            if exchange["type"] == "technosphere":
                exc_id = (
                    exchange["name"],
                    exchange.get("reference product"),
                    exchange["unit"],
                )

                if exc_id in tuples_of_blacklisted_datasets:
                    for d in blacklist:
                        if exc_id == (d["name"], d.get("reference product"), d["unit"]):
                            if "replacement" in d:
                                exchange["name"] = d["replacement"]["name"]
                                exchange["reference product"] = d["replacement"][
                                    "reference product"
                                ]
                                exchange["location"] = d["replacement"]["location"]

    return datasets


def check_amount_format(database: list) -> list:
    """
    Check that the `amount` field is of type `float`.
    :param database: database to check
    :return: database with corrected amount field
    """

    for dataset in database:
        for exc in dataset["exchanges"]:
            if not isinstance(exc["amount"], float):
                exc["amount"] = float(exc["amount"])

            if isinstance(exc["amount"], (np.float64, np.ndarray)):
                exc["amount"] = float(exc["amount"])

        for k, v in dataset.items():
            if isinstance(v, dict):
                for i, j in v.items():
                    if isinstance(j, (np.float64, np.ndarray)):
                        v[i] = float(v[i])

        for e in dataset["exchanges"]:
            for k, v in e.items():
                if isinstance(v, (np.float64, np.ndarray)):
                    e[k] = float(e[k])

    return database


class BaseInventoryImport:
    """
    Base class for inventories that are to be merged with the wurst database.

    :ivar database: the target database for the import (the ecoinvent database),
    unpacked to a list of dicts
    :ivar version_in: the ecoinvent database version of the inventory to import
    :ivar version_out: the ecoinvent database version the imported inventories
    should comply with
    :ivar path: the filepath of the inventories to import

    """

    def __init__(
        self,
        database: List[dict],
        version_in: str,
        version_out: str,
        path: Union[str, Path],
        system_model: str,
        keep_uncertainty_data: bool = False,
    ) -> None:
        """Create a :class:`BaseInventoryImport` instance."""
        self.database = database
        self.db_code = [x["code"] for x in self.database]
        self.db_names = [
            (x["name"], x["reference product"], x["location"]) for x in self.database
        ]
        self.version_in = version_in
        self.version_out = version_out
        self.biosphere_dict = get_biosphere_code(self.version_out)
        self.correspondence_bio_flows = get_correspondence_bio_flows()
        self.system_model = system_model
        self.consequential_blacklist = get_consequential_blacklist()
        self.list_unlinked = []
        self.keep_uncertainty_data = keep_uncertainty_data

        if "http" in str(path):
            r = requests.head(path)
            if r.status_code != 200:
                raise ValueError("The file at {} could not be found.".format(path))
        else:
            if not Path(path).exists():
                raise FileNotFoundError(
                    f"The inventory file {path} could not be found."
                )

        self.path = Path(path) if isinstance(path, str) else path
        self.import_db = self.load_inventory(path)

        # register migration maps
        # as imported inventories link
        # to different ecoinvent versions
        ei_versions = ["35", "36", "37", "38", "39"]

        for combination in itertools.product(ei_versions, ei_versions):
            if combination[0] != combination[1]:
                mapping = generate_migration_maps(combination[0], combination[1])
                if len(mapping["data"]) > 0:
                    Migration(f"migration_{combination[0]}_{combination[1]}").write(
                        mapping,
                        description=f"Change technosphere names due to change from {combination[0]} to {combination[1]}",
                    )

    def load_inventory(self, path: Union[str, Path]) -> None:
        """Load an inventory from a specified path.
        Sets the :attr:`import_db` attribute.
        :param str path: Path to the inventory file
        :returns: Nothing.
        """
        return None

    def prepare_inventory(self) -> None:
        """Prepare the inventory for the merger with Ecoinvent.
        Modifies :attr:`import_db` in-place.
        :returns: Nothing
        """

    def check_for_already_existing_datasets(self) -> None:
        """
        Check whether the inventories to be imported are not
        already in the source database.
        """

        # print if we find datasets that already exist
        already_exist = [
            (x["name"].lower(), x["reference product"].lower(), x["location"])
            for x in self.import_db.data
            if x["code"] in self.db_code
        ]

        already_exist.extend(
            [
                (x["name"].lower(), x["reference product"].lower(), x["location"])
                for x in self.import_db.data
                if (x["name"].lower(), x["reference product"].lower(), x["location"])
                in self.db_names
            ]
        )

        if len(already_exist) > 0:
            print(
                "The following datasets to import already exist "
                "in the source database. "
                "They will not be imported"
            )
            table = PrettyTable(["Name", "Reference product", "Location", "File"])

            if isinstance(self.path, str):
                name = self.path
            else:
                name = self.path.name

            for dataset in already_exist:
                table.add_row([dataset[0][:50], dataset[1][:30], dataset[2], name])

            print(table)

        self.import_db.data = [
            x for x in self.import_db.data if x["code"] not in self.db_code
        ]
        self.import_db.data = [
            x
            for x in self.import_db.data
            if (x["name"], x["reference product"], x["location"]) not in self.db_names
        ]

    def merge_inventory(self) -> List[dict]:
        """Prepare :attr:`import_db` and merge the inventory to the ecoinvent :attr:`database`.
        Calls :meth:`prepare_inventory`. Changes the :attr:`database` attribute.
        :returns: Nothing
        """

        self.prepare_inventory()
        return self.import_db

    def search_missing_exchanges(self, label: str, value: str) -> List[dict]:
        """
        Return a list of activities for which
        a given exchange cannot be found
        :param label: the label of the field to look for
        :param value: the value of the field to look for
        :return:
        """

        results = []
        for act in self.import_db.data:
            if (
                len([a for a in act["exchanges"] if label in a and a[label] == value])
                == 0
            ):
                results.append(act)

        return results

    def search_missing_field(self, field: str, scope: str = "activity") -> List[dict]:
        """Find exchanges and activities that do not contain a specific field
        in :attr:`imort_db`
        :param str field: label of the field to search for.
        :param scope: "activity" or "all". whether to search in the activity
        or the activity and its exchanges
        :returns: a list of dictionaries, activities and exchanges
        :rtype: list
        """
        results = []
        for act in self.import_db.data:
            if field not in act:
                results.append(act)

            if scope == "all":
                for ex in act["exchanges"]:
                    if ex["type"] == "technosphere" and field not in ex:
                        results.append(ex)
        return results

    def add_product_field_to_exchanges(self) -> None:
        """Add the `product` key to the production and
        technosphere exchanges in :attr:`import_db`.
        Also add `code` field if missing.
        For production exchanges, use the value of the `reference_product` field.
        For technosphere exchanges, search the activities in :attr:`import_db` and
        use the reference product. If none is found, search the Ecoinvent :attr:`database`.
        Modifies the :attr:`import_db` attribute in place.
        :raises IndexError: if no corresponding activity (and reference product) can be found.
        """
        # Add a `product` field to the production exchange
        for dataset in self.import_db.data:
            for exchange in dataset["exchanges"]:
                if exchange["type"] == "production":
                    if "product" not in exchange:
                        exchange["product"] = dataset["reference product"]

                    if exchange["name"] != dataset["name"]:
                        exchange["name"] = dataset["name"]

        # Add a `product` field to technosphere exchanges
        for dataset in self.import_db.data:
            for exchange in dataset["exchanges"]:
                if exchange["type"] == "technosphere":
                    # Check if the field 'product' is present
                    if not "product" in exchange:
                        exchange["product"] = self.correct_product_field(
                            (
                                exchange["name"],
                                exchange["location"],
                                exchange["unit"],
                                exchange.get("reference product", None),
                            )
                        )

                    # If a 'reference product' field is present, we make sure
                    # it matches with the new 'product' field
                    # if "reference product" in y:
                    if "reference product" in exchange:
                        try:
                            assert exchange["product"] == exchange["reference product"]
                        except AssertionError:
                            exchange["product"] = self.correct_product_field(
                                (
                                    exchange["name"],
                                    exchange["location"],
                                    exchange["unit"],
                                    exchange.get("reference product", None),
                                )
                            )

        # Add a `code` field if missing
        for dataset in self.import_db.data:
            if "code" not in dataset:
                dataset["code"] = str(uuid.uuid4().hex)

    @lru_cache
    def correct_product_field(self, exc: tuple) -> [str, None]:
        """
        Find the correct name for the `product` field of the exchange
        :param exc: a dataset exchange
        :return: name of the product field of the exchange

        """
        # Look first in the imported inventories
        candidate = next(
            ws.get_many(
                self.import_db.data,
                ws.equals("name", exc[0]),
                ws.equals("location", exc[1]),
                ws.equals("unit", exc[2]),
            ),
            None,
        )

        # If not, look in the ecoinvent inventories
        if candidate is None:
            if exc[-1] is not None:
                candidate = next(
                    ws.get_many(
                        self.database,
                        ws.equals("name", exc[0]),
                        ws.equals("location", exc[1]),
                        ws.equals("unit", exc[2]),
                        ws.equals("reference product", exc[-1]),
                    ),
                    None,
                )
            else:
                candidate = next(
                    ws.get_many(
                        self.database,
                        ws.equals("name", exc[0]),
                        ws.equals("location", exc[1]),
                        ws.equals("unit", exc[2]),
                    ),
                    None,
                )

        if candidate is not None:
            return candidate["reference product"]

        self.list_unlinked.append(
            (
                exc[0],
                exc[-1],
                exc[1],
                None,
                exc[2],
                "technosphere",
                self.path.name,
            )
        )

        return None

    def add_biosphere_links(self, delete_missing: bool = False) -> None:
        """Add links for biosphere exchanges to :attr:`import_db`
        Modifies the :attr:`import_db` attribute in place.

        :param delete_missing: whether unlinked exchanges should be deleted or not.
        """
        for x in self.import_db.data:
            for y in x["exchanges"]:
                if y["type"] == "biosphere":
                    if isinstance(y["categories"], str):
                        y["categories"] = tuple(y["categories"].split("::"))

                    if len(y["categories"]) > 1:
                        key = (
                            y["name"],
                            y["categories"][0],
                            y["categories"][1],
                            y["unit"],
                        )
                    else:
                        key = (
                            y["name"],
                            y["categories"][0].strip(),
                            "unspecified",
                            y["unit"],
                        )

                    if key not in self.biosphere_dict:
                        if self.correspondence_bio_flows.get(key[1], {}).get(key[0]):
                            new_key = list(key)
                            new_key[0] = self.correspondence_bio_flows[key[1]][key[0]]
                            key = tuple(new_key)
                            assert (
                                key in self.biosphere_dict
                            ), f"Could not find a biosphere flow for {key}."
                        else:
                            print(key)
                            continue

                    y["input"] = (
                        "biosphere3",
                        self.biosphere_dict[key],
                    )

    def lower_case_technosphere_exchanges(self) -> None:
        blakclist = [
            "NOx",
            "SOx",
            "N-",
            "EUR",
        ]

        for ds in self.import_db.data:
            # lower case name and reference product
            if not any([x in ds["name"] for x in blakclist]):
                ds["name"] = ds["name"][0].lower() + ds["name"][1:]
            if not any([x in ds["reference product"] for x in blakclist]):
                ds["reference product"] = (
                    ds["reference product"][0].lower() + ds["reference product"][1:]
                )

            for exc in ds["exchanges"]:
                if exc["type"] in ["technosphere", "production"]:
                    if not any([x in exc["name"] for x in blakclist]):
                        exc["name"] = exc["name"][0].lower() + exc["name"][1:]

                    if not any(
                        [x in exc.get("reference product", "") for x in blakclist]
                    ):
                        if exc.get("reference product") is not None:
                            exc["reference product"] = (
                                exc["reference product"][0].lower()
                                + exc["reference product"][1:]
                            )

                    if not any([x in exc.get("product", "") for x in blakclist]):
                        if exc.get("product") is not None:
                            exc["product"] = (
                                exc["product"][0].lower() + exc["product"][1:]
                            )

    def remove_ds_and_modifiy_exchanges(self, name: str, ex_data: dict) -> None:
        """
        Remove an activity dataset from :attr:`import_db` and replace the corresponding
        technosphere exchanges by what is given as second argument.
        :param str name: name of activity to be removed
        :param dict ex_data: data to replace the corresponding exchanges
        :returns: Nothing
        """

        self.import_db.data = [
            act for act in self.import_db.data if not act["name"] == name
        ]

        for act in self.import_db.data:
            for ex in act["exchanges"]:
                if ex["type"] == "technosphere" and ex["name"] == name:
                    ex.update(ex_data)
                    # make sure there is no existing link
                    if "input" in ex:
                        del ex["input"]

            # Delete any field that does not have information
            for key in act:
                if act[key] is None:
                    act.pop(key)

    def display_unlinked_exchanges(self):
        """
        Display the list of unlinked exchanges
        using prettytable
        """
        print("List of unlinked exchanges:")

        table = PrettyTable()
        table.field_names = [
            "Name",
            "Reference product",
            "Location",
            "Categories",
            "Unit",
            "Type",
            "File",
        ]
        table.add_rows(list(set(self.list_unlinked)))
        print(table)


class DefaultInventory(BaseInventoryImport):
    """
    Importing class. Inherits from :class:`BaseInventoryImport`.

    """

    def __init__(
        self,
        database,
        version_in,
        version_out,
        path,
        system_model,
        keep_uncertainty_data,
    ):
        super().__init__(
            database, version_in, version_out, path, system_model, keep_uncertainty_data
        )

    def load_inventory(self, path: Union[str, Path]) -> bw2io.ExcelImporter:
        return ExcelImporter(path)

    def prepare_inventory(self) -> None:
        if self.version_in != self.version_out:
            # if version_out is 3.9, migrate towards 3.8 first, then 3.9
            if self.version_out in ["3.9", "3.9.1"]:
                print("Migrating to 3.8 first")
                if self.version_in != "3.8":
                    self.import_db.migrate(
                        f"migration_{self.version_in.replace('.', '')}_38"
                    )
                self.import_db.migrate(
                    f"migration_38_{self.version_out.replace('.', '')}"
                )
            self.import_db.migrate(
                f"migration_{self.version_in.replace('.', '')}_{self.version_out.replace('.', '')}"
            )

        if self.system_model == "consequential":
            self.import_db.data = (
                check_for_datasets_compliance_with_consequential_database(
                    self.import_db.data, self.consequential_blacklist
                )
            )

        self.import_db.data = remove_categories(self.import_db.data)

        self.lower_case_technosphere_exchanges()
        self.add_biosphere_links()
        self.add_product_field_to_exchanges()

        # Remove uncertainty data
        if not self.keep_uncertainty_data:
            print("Remove uncertainty data.")
            self.database = remove_uncertainty(self.database)

        # Check for duplicates
        self.check_for_already_existing_datasets()
        self.import_db.data = check_for_duplicate_datasets(self.import_db.data)

        if self.list_unlinked:
            self.display_unlinked_exchanges()


class VariousVehicles(BaseInventoryImport):
    """
    Imports various future vehicles' inventories (two-wheelers, buses, trams, etc.).

    :ivar database: wurst database
    :ivar version_in: original ecoinvent version of the inventories
    :ivar version_out: ecoinvent version the inventories should comply with
    :ivar path: filepath of the inventories
    :ivar year: year of the database
    :ivar regions:
    :ivar model: IAM model
    :ivar scenario: IAM scenario
    :ivar vehicle_type: "two-wheeler", "car, "truck" or "bus"
    :ivar relink: whether suppliers within a dataset need to be relinked
    :ivar has_fleet: whether the `vehicle_type` has associated fleet information
    """

    def __init__(
        self,
        database: List[dict],
        version_in: str,
        version_out: str,
        path: Union[str, Path],
        year: int,
        regions: List[str],
        model: str,
        scenario: str,
        vehicle_type: str,
        relink: bool = False,
        has_fleet: bool = False,
        system_model: str = "cutoff",
    ) -> None:
        super().__init__(database, version_in, version_out, path, system_model)
        self.year = year
        self.regions = regions
        self.model = model
        self.scenario = scenario
        self.vehicle_type = vehicle_type
        self.relink = relink
        self.has_fleet = has_fleet
        self.geo = Geomap(model=model)

    def load_inventory(self, path):
        return ExcelImporter(path)

    def prepare_inventory(self):
        # if version_out is 3.9, migrate towards 3.8 first, then 3.9
        if self.version_out in ["3.9", "3.9.1"]:
            print("Migrating to 3.8 first")
            if self.version_in != "3.8":
                self.import_db.migrate(
                    f"migration_{self.version_in.replace('.', '')}_38"
                )
            self.import_db.migrate(f"migration_38_{self.version_out.replace('.', '')}")
        self.import_db.migrate(
            f"migration_{self.version_in.replace('.', '')}_{self.version_out.replace('.', '')}"
        )

        self.lower_case_technosphere_exchanges()
        self.add_biosphere_links()
        self.add_product_field_to_exchanges()
        # Check for duplicates
        self.check_for_already_existing_datasets()

        if self.list_unlinked:
            self.display_unlinked_exchanges()

    def merge_inventory(self):
        self.database.extend(self.import_db.data)

        print("Done!")

        return self.database


class AdditionalInventory(BaseInventoryImport):
    """
    Import additional inventories, if any.
    """

    def __init__(self, database, version_in, version_out, path, system_model):
        super().__init__(database, version_in, version_out, path, system_model)

    def load_inventory(self, path):
        if "http" in path:
            # online file
            # we need to save it locally first
            response = requests.get(path)
            Path(DATA_DIR / "cache").mkdir(parents=True, exist_ok=True)
            path = str(Path(DATA_DIR / "cache" / "temp.csv"))
            with open(path, "w", encoding="utf-8") as f:
                writer = csv.writer(
                    f,
                    quoting=csv.QUOTE_NONE,
                    delimiter=",",
                    quotechar="'",
                    escapechar="\\",
                )
                for line in response.iter_lines():
                    writer.writerow(line.decode("utf-8").split(","))

        if Path(path).suffix == ".xlsx":
            return ExcelImporter(path)
        elif Path(path).suffix == ".csv":
            return CSVImporter(path)
        else:
            raise ValueError(
                "Incorrect filetype for inventories." "Should be either .xlsx or .csv"
            )

    def remove_missing_fields(self):
        """
        Remove any field that does not have information.
        """

        for dataset in self.import_db.data:
            for key, value in list(dataset.items()):
                if not value:
                    del dataset[key]

    def prepare_inventory(self):
        if str(self.version_in) != self.version_out:
            # if version_out is 3.9, migrate towards 3.8 first, then 3.9
            if self.version_out in ["3.9", "3.9.1"]:
                if str(self.version_in) != "3.8":
                    print("Migrating to 3.8 first")
                    self.import_db.migrate(
                        f"migration_{self.version_in.replace('.', '')}_38"
                    )
                self.import_db.migrate(
                    f"migration_38_{self.version_out.replace('.', '')}"
                )

            self.import_db.migrate(
                f"migration_{self.version_in.replace('.', '')}_{self.version_out.replace('.', '')}"
            )

        if self.system_model == "consequential":
            self.import_db.data = (
                check_for_datasets_compliance_with_consequential_database(
                    self.import_db.data, self.consequential_blacklist
                )
            )

        list_missing_prod = self.search_missing_exchanges(
            label="type", value="production"
        )

        if len(list_missing_prod) > 0:
            print("The following datasets are missing a `production` exchange.")
            print("You should fix those before proceeding further.\n")
            table = PrettyTable(
                ["Name", "Reference product", "Location", "Unit", "File"]
            )
            for dataset in list_missing_prod:
                table.add_row(
                    [
                        dataset.get("name", "XXXX"),
                        dataset.get("referece product", "XXXX"),
                        dataset.get("location", "XXXX"),
                        dataset.get("unit", "XXXX"),
                        self.path.name,
                    ]
                )

            print(table)

            sys.exit()

        self.import_db.data = remove_categories(self.import_db.data)
        self.add_biosphere_links(delete_missing=True)
        list_missing_ref = self.search_missing_field(field="name")
        list_missing_ref.extend(self.search_missing_field(field="reference product"))
        list_missing_ref.extend(self.search_missing_field(field="location"))
        list_missing_ref.extend(self.search_missing_field(field="unit"))

        if len(list_missing_ref) > 0:
            print(
                "The following datasets are missing an important field "
                "(`name`, `reference product`, `location` or `unit`).\n"
            )

            print("You should fix those before proceeding further.\n")
            table = PrettyTable(
                ["Name", "Reference product", "Location", "Unit", "File"]
            )
            for dataset in list_missing_ref:
                table.add_row(
                    [
                        dataset.get("name", "XXXX"),
                        dataset.get("referece product", "XXXX"),
                        dataset.get("location", "XXXX"),
                        dataset.get("unit", "XXXX"),
                        self.path.name,
                    ]
                )

            print(table)

        if len(list_missing_prod) > 0 or len(list_missing_ref) > 0:
            sys.exit()

        self.remove_missing_fields()
        self.add_product_field_to_exchanges()
        # Check for duplicates
        self.check_for_already_existing_datasets()
        self.import_db.data = check_for_duplicate_datasets(self.import_db.data)
        # check numbers format
        self.import_db.data = check_amount_format(self.import_db.data)

        if self.list_unlinked:
            self.display_unlinked_exchanges()
