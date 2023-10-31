import itertools
from copy import copy

from bw2data import Database, databases
from bw2io.importers.base_lci import LCIImporter
from wurst.linking import change_db_name, check_internal_linking, link_internal


class BW25Importer(LCIImporter):
    def __init__(self, db_name, data):
        self.db_name = db_name
        self.data = data
        for act in self.data:
            act["database"] = self.db_name

    # we override `write_database`
    # to allow existing databases
    # to be overwritten
    def write_database(self):
        def no_exchange_generator(data):
            for ds in data:
                cp = copy(ds)
                cp["exchanges"] = []
                yield cp

        if self.db_name in databases:
            print(f"Database {self.db_name} already exists: " "it will be overwritten.")
        super().write_database(
            list(no_exchange_generator(self.data)), backend="iotable"
        )

        dependents = {exc["input"][0] for ds in self.data for exc in ds["exchanges"]}
        lookup = {
            obj.key: obj.id
            for obj in itertools.chain(*[Database(label) for label in dependents])
        }

        def technosphere_generator(data, lookup):
            for ds in data:
                target = lookup[(ds["database"], ds["code"])]
                for exc in ds["exchanges"]:
                    if exc["type"] in (
                        "substitution",
                        "production",
                        "generic production",
                    ):
                        yield {
                            "row": lookup[exc["input"]],
                            "col": target,
                            "amount": exc["amount"],
                            "flip": False,
                        }
                    elif exc["type"] == "technosphere":
                        yield {
                            "row": lookup[exc["input"]],
                            "col": target,
                            "amount": exc["amount"],
                            "flip": True,
                        }

        def biosphere_generator(data, lookup):
            for ds in data:
                target = lookup[(ds["database"], ds["code"])]
                for exc in ds["exchanges"]:
                    if exc["type"] == "biosphere":
                        yield {
                            "row": lookup[exc["input"]],
                            "col": target,
                            "amount": exc["amount"],
                            "flip": False,
                        }

        Database(self.db_name).write_exchanges(
            technosphere_generator(self.data, lookup),
            biosphere_generator(self.data, lookup),
            list(dependents),
        )


def write_brightway_database(data, name):
    # Restore parameters to Brightway2 format
    # which allows for uncertainty and comments
    change_db_name(data, name)
    link_internal(data)
    check_internal_linking(data)
    BW25Importer(name, data).write_database()
