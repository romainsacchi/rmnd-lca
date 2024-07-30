"""
This module contains functions to write a Brightway 2.5 database.
"""

from bw2data import Database, databases
from bw2io.importers.base_lci import LCIImporter
from wurst.linking import change_db_name, check_internal_linking, link_internal


class BW25Importer(LCIImporter):
    """
    Class to write a Brightway 2.5 database from a Wurst database.
    """

    def __init__(self, db_name: str, data: list) -> None:
        """
        :param db_name: Name of the database to write
        :type db_name: str
        :param data: Wurst database
        :type data: list
        """
        super().__init__(db_name)
        self.db_name = db_name
        self.data = data
        for act in self.data:
            act["database"] = self.db_name


def write_brightway_database(data: list, name: str) -> None:
    # Restore parameters to Brightway2 format
    # which allows for uncertainty and comments
    change_db_name(data, name)
    link_internal(data)
    check_internal_linking(data)

    if name in databases:
        print(f"Database {name} already exists: it will be overwritten.")
        del databases[name]

    BW25Importer(name, data).write_database()
