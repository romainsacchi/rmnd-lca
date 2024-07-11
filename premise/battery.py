"""
module to adjust the battery inputs to reflect progress in
terms of cell energy density.

"""

import yaml

from .filesystem_constants import DATA_DIR
from .logger import create_logger
from .transformation import BaseTransformation, IAMDataCollection, List, np, ws

logger = create_logger("battery")


def load_cell_energy_density():
    """
    Load cell energy density data.
    """
    with open(DATA_DIR / "battery/energy_density.yaml", "r") as file:
        data = yaml.load(file, Loader=yaml.FullLoader)

    return {x["ecoinvent_aliases"]["name"]: x["target"] for x in data.values()}


def _update_battery(scenario, version, system_model):
    battery = Battery(
        database=scenario["database"],
        iam_data=scenario["iam data"],
        model=scenario["model"],
        pathway=scenario["pathway"],
        year=scenario["year"],
        version=version,
        system_model=system_model,
        cache=scenario.get("cache"),
        index=scenario.get("index"),
    )

    battery.adjust_battery_mass()

    scenario["database"] = battery.database
    scenario["index"] = battery.index
    scenario["cache"] = battery.cache

    return scenario


class Battery(BaseTransformation):
    """
    Class that modifies the battery market to reflect progress
    in terms of cell energy density.

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
        super().__init__(
            database,
            iam_data,
            model,
            pathway,
            year,
            version,
            system_model,
            cache,
            index,
        )
        self.system_model = system_model

    def adjust_battery_mass(self) -> None:
        """
        Adjust vehicle components (e.g., battery).
        Adjust the battery mass to reflect progress in battery technology.
        Specifically, we adjust the battery mass to reflect progress in
        terms of cell energy density.
        We leave the density unchanged after 2050.
        """

        energy_density = load_cell_energy_density()

        filters = [ws.contains("name", x) for x in energy_density]

        for ds in ws.get_many(
            self.database,
            ws.exclude(
                ws.either(
                    *[
                        ws.contains("name", x)
                        for x in [
                            "market for battery",
                            "battery production",
                            "battery cell production",
                            "cell module production",
                        ]
                    ]
                )
            ),
        ):

            for exc in ws.technosphere(ds, ws.either(*filters)):
                name = [x for x in energy_density if x in exc["name"]][0]

                scaling_factor = energy_density[name][2020] / np.clip(
                    np.interp(
                        self.year,
                        list(energy_density[name].keys()),
                        list(energy_density[name].values()),
                    ),
                    0.1,
                    0.5,
                )

                if "log parameters" not in ds:
                    ds["log parameters"] = {}

                ds["log parameters"]["battery input"] = exc["name"]
                ds["log parameters"]["old battery mass"] = exc["amount"]
                exc["amount"] *= scaling_factor
                ds["log parameters"]["new battery mass"] = exc["amount"]

                self.write_log(ds, status="modified")

    def write_log(self, dataset, status="created"):
        """
        Write log file.
        """

        logger.info(
            f"{status}|{self.model}|{self.scenario}|{self.year}|"
            f"{dataset['name']}|{dataset['location']}|"
            f"{dataset.get('log parameters', {}).get('battery input', '')}|"
            f"{dataset.get('log parameters', {}).get('old battery mass', '')}|"
            f"{dataset.get('log parameters', {}).get('new battery mass', '')}"
        )
