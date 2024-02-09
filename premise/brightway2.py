"""
Class to write a Brightway2 database from a Wurst database.
"""

from bw2data import databases
from bw2io.importers.base_lci import LCIImporter
from wurst.linking import change_db_name, check_internal_linking, link_internal


class BW2Importer(LCIImporter):
    """
    Class to write a Brightway2 database from a Wurst database.
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

    # we override `write_database`
    # to allow existing databases
    # to be overwritten
    def write_database(self):
        if self.db_name in databases:
            print(f"Database {self.db_name} already exists: it will be overwritten.")
            del databases[self.db_name]
        super().write_database()


def write_brightway_database(data: list, name: str) -> None:
    """
    Write a Brightway2 database from a Wurst database.
    """
    # Restore parameters to Brightway2 format
    # which allows for uncertainty and comments
    change_db_name(data, name)
    link_internal(data)
    check_internal_linking(data)
    BW2Importer(name, data).write_database()
