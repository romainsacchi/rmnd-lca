"""
module to adjust the battery inputs to reflect progress in
terms of cell energy density.

"""

import yaml

from .filesystem_constants import DATA_DIR
from .logger import create_logger
from .transformation import BaseTransformation, IAMDataCollection, List, np, ws
from .utils import rescale_exchanges

logger = create_logger("battery")


def load_cell_energy_density():
    """
    Load cell energy density data.
    """
    with open(DATA_DIR / "battery/energy_density.yaml", "r") as file:
        data = yaml.load(file, Loader=yaml.FullLoader)

    result = {}
    for key, value in data.items():
        names = value["ecoinvent_aliases"]["name"]
        if isinstance(names, list):
            for name in names:
                result[name] = value["target"]
        else:
            result[names] = value["target"]

    return result


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

    if battery.iam_data.battery_scenarios is not None:
        battery.adjust_battery_market_shares()

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

    def adjust_battery_market_shares(self) -> None:
        """
        Based on scenario data, adjust the shares within the datasets:
        - market for battery capacity (MIX scenario)
        - market for battery capacity (LFP scenario)
        - market for battery capacity (NCx scenario)
        - market for battery capacity (PLiB scenario)
        """

        market_datasets = {
            "market for battery capacity (MIX scenario)": "MIX",
            "market for battery capacity (LFP scenario)": "LFP",
            "market for battery capacity (NCx scenario)": "NCX",
            "market for battery capacity (PLiB scenario)": "PLIB",
        }

        datasets_mapping = {
            v: k
            for k, v in {
                "LAB": "market for battery capacity, Li-ion, Li-O2",
                "LFP": "market for battery capacity, Li-ion, LFP",
                "LSB": "market for battery capacity, Li-sulfur, Li-S",
                "NCA": "market for battery capacity, Li-ion, NCA",
                "NMC111": "market for battery capacity, Li-ion, NMC111",
                "NMC532": "market for battery capacity, Li-ion, NMC523",
                "NMC622": "market for battery capacity, Li-ion, NMC622",
                "NMC811": "market for battery capacity, Li-ion, NMC811",
                "NMC900-Si": "market for battery capacity, Li-ion, NMC955",
                "SIB": "market for battery capacity, Sodium-ion, SiB",
            }.items()
        }

        for ds in ws.get_many(
            self.database,
            ws.either(*[ws.equals("name", name) for name in market_datasets]),
        ):

            if self.year in self.iam_data.battery_scenarios.year:
                shares = self.iam_data.battery_scenarios.sel(
                    scenario=market_datasets[ds["name"]],
                    year=self.year,
                )
            elif self.year < min(self.iam_data.battery_scenarios.year):
                shares = self.iam_data.battery_scenarios.sel(
                    scenario=market_datasets[ds["name"]],
                    year=min(self.iam_data.battery_scenarios.year),
                )
            elif self.year > max(self.iam_data.battery_scenarios.year):
                shares = self.iam_data.battery_scenarios.sel(
                    scenario=market_datasets[ds["name"]],
                    year=max(self.iam_data.battery_scenarios.year),
                )
            else:
                shares = self.iam_data.battery_scenarios.sel(
                    scenario=market_datasets[ds["name"]],
                ).interp(year=self.year)

            if "log parameters" not in ds:
                ds["log parameters"] = {}

            for exc in ws.technosphere(ds):
                if exc["name"] in datasets_mapping:
                    exc["amount"] = shares.sel(
                        chemistry=datasets_mapping[exc["name"]]
                    ).values.item()

                    ds["log parameters"][
                        f"{datasets_mapping[exc['name']]} market share"
                    ] = exc["amount"]

            self.write_log(ds, status="modified")

    def adjust_battery_mass(self) -> None:
        """
        Adjust vehicle components (e.g., battery).
        Adjust the battery mass to reflect progress in battery technology.
        Specifically, we adjust the battery mass to reflect progress in
        terms of cell energy density.
        We leave the density unchanged after 2050.
        """

        energy_density = load_cell_energy_density()

        for ds in ws.get_many(
            self.database,
            ws.contains("name", "market for battery capacity"),
        ):
            if ds["name"] in energy_density:

                mean_2020_energy_density = energy_density[ds["name"]][2020]["mean"]
                minimum_2020_energy_density = energy_density[ds["name"]][2020][
                    "minimum"
                ]
                maximum_2020_energy_density = energy_density[ds["name"]][2020][
                    "maximum"
                ]
                mean_2050_energy_density = energy_density[ds["name"]][2050]["mean"]
                minimum_2050_energy_density = energy_density[ds["name"]][2050][
                    "minimum"
                ]
                maximum_2050_energy_density = energy_density[ds["name"]][2050][
                    "maximum"
                ]

                scaling_factor = mean_2020_energy_density / np.clip(
                    np.interp(
                        self.year,
                        [2020, 2050],
                        [mean_2020_energy_density, mean_2050_energy_density],
                    ),
                    0,
                    1,
                )

                scaling_factor_min = minimum_2020_energy_density / np.clip(
                    np.interp(
                        self.year,
                        [2020, 2050],
                        [minimum_2020_energy_density, minimum_2050_energy_density],
                    ),
                    0,
                    1,
                )

                scaling_factor_max = maximum_2020_energy_density / np.clip(
                    np.interp(
                        self.year,
                        [2020, 2050],
                        [maximum_2020_energy_density, maximum_2050_energy_density],
                    ),
                    0,
                    1,
                )

                if "log parameters" not in ds:
                    ds["log parameters"] = {}

                ds["log parameters"]["battery input"] = [
                    e["name"]
                    for e in ws.technosphere(
                        ds, ws.contains("name", "market for battery")
                    )
                ][0]

                ds["log parameters"]["old battery mass"] = sum(
                    e["amount"]
                    for e in ws.technosphere(
                        ds, ws.contains("name", "market for battery")
                    )
                )

                for exc in ws.technosphere(ds, ws.equals("unit", "kilogram")):
                    exc["amount"] *= scaling_factor
                    exc["minimum"] *= scaling_factor_min
                    exc["maximum"] *= scaling_factor_max

                ds["log parameters"]["new battery mass"] = sum(
                    e["amount"]
                    for e in ws.technosphere(
                        ds, ws.contains("name", "market for battery")
                    )
                )

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
            f"{dataset.get('log parameters', {}).get('new battery mass', '')}|"
            f"{dataset.get('log parameters', {}).get('NMC111 market share', '')}|"
            f"{dataset.get('log parameters', {}).get('NMC532 market share', '')}|"
            f"{dataset.get('log parameters', {}).get('NMC622 market share', '')}|"
            f"{dataset.get('log parameters', {}).get('NMC811 market share', '')}|"
            f"{dataset.get('log parameters', {}).get('NMC900-Si market share', '')}|"
            f"{dataset.get('log parameters', {}).get('LFP market share', '')}|"
            f"{dataset.get('log parameters', {}).get('NCA market share', '')}|"
            f"{dataset.get('log parameters', {}).get('LAB market share', '')}|"
            f"{dataset.get('log parameters', {}).get('LSB market share', '')}|"
            f"{dataset.get('log parameters', {}).get('SIB market share', '')}"
        )
