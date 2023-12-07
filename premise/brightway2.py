from bw2data import databases
from bw2io.importers.base_lci import LCIImporter
from wurst.linking import change_db_name, check_internal_linking, link_internal

from .utils import reset_all_codes


class BW2Importer(LCIImporter):
    def __init__(self, db_name, data):
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


def write_brightway_database(data, name):
    # Restore parameters to Brightway2 format
    # which allows for uncertainty and comments
    change_db_name(data, name)
    link_internal(data)
    check_internal_linking(data)
    BW2Importer(name, data).write_database()
