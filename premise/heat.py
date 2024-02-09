"""
Integrates projections regarding heat production and supply.
"""

from .logger import create_logger
from .transformation import BaseTransformation, IAMDataCollection, List, ws

logger = create_logger("heat")


def _update_heat(scenario, version, system_model):
    heat = Heat(
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

    heat.fetch_fuel_market_co2_emissions()
    heat.regionalize_heat_production()
    heat.relink_datasets()
    scenario["database"] = heat.database
    scenario["cache"] = heat.cache
    scenario["index"] = heat.index

    return scenario


class Heat(BaseTransformation):
    """
    Class that modifies fuel inventories and markets
    in ecoinvent based on IAM output data.
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
    ):
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

        self.carbon_intensity_markets = {}

    def fetch_fuel_market_co2_emissions(self):
        """
        Fetch CO2 emissions from fuel markets.
        """
        fuel_markets = [
            "market for diesel, low-sulfur",
            "market for petrol, low-sulfur",
            "market for natural gas, high pressure",
        ]

        for dataset in ws.get_many(
            self.database,
            ws.either(*[ws.equals("name", n) for n in fuel_markets]),
        ):
            if "log parameters" in dataset:
                if "fossil CO2 per kg fuel" in dataset["log parameters"]:
                    self.carbon_intensity_markets[
                        (dataset["name"], dataset["location"])
                    ] = {"fossil": dataset["log parameters"]["fossil CO2 per kg fuel"]}
                if "non-fossil CO2 per kg fuel" in dataset["log parameters"]:
                    self.carbon_intensity_markets[
                        (dataset["name"], dataset["location"])
                    ].update(
                        {
                            "non-fossil": dataset["log parameters"][
                                "non-fossil CO2 per kg fuel"
                            ]
                        }
                    )

        # add "market for diesel" to self.carbon_intensity_markets
        # by duplicating the "market for low-sulfur" entries
        # add "market for natural gas, low pressure" to self.carbon_intensity_markets
        # by duplicating the "market for natural gas, high pressure" entries

        new_keys = {}
        for key, value in self.carbon_intensity_markets.items():
            if key[0] == "market for diesel, low-sulfur":
                new_keys[("market for diesel", key[1])] = value
                new_keys[("market group for diesel", key[1])] = value
                new_keys[("market group for diesel, low-sulfur", key[1])] = value
            if key[0] == "market for petrol, low-sulfur":
                new_keys[("market for petrol", key[1])] = value
                new_keys[("market for petrol, unleaded", key[1])] = value
            if key[0] == "market for natural gas, high pressure":
                new_keys[("market group for natural gas, high pressure", key[1])] = (
                    value
                )
                new_keys[("market for natural gas, low pressure", key[1])] = value

        self.carbon_intensity_markets.update(new_keys)

    def regionalize_heat_production(self):
        """
        Regionalize heat production.

        """

        for heat_datasets in self.heat_techs.values():
            for dataset in ws.get_many(
                self.database,
                ws.either(*[ws.contains("name", n) for n in heat_datasets]),
                ws.equals("unit", "megajoule"),
                ws.doesnt_contain_any("location", self.regions),
            ):
                new_ds = self.fetch_proxies(
                    name=dataset["name"],
                    ref_prod=dataset["reference product"],
                    exact_name_match=True,
                    exact_product_match=True,
                )

                if len(new_ds) == 0:
                    continue

                for ds in new_ds.values():
                    fossil_co2, non_fossil_co2 = 0.0, 0.0

                    for exc in ws.technosphere(ds):
                        if (
                            exc["name"],
                            exc["location"],
                        ) in self.carbon_intensity_markets:
                            fossil_co2 += exc[
                                "amount"
                            ] * self.carbon_intensity_markets.get(
                                (exc["name"], exc["location"]), {}
                            ).get(
                                "fossil", 0.0
                            )
                            non_fossil_co2 += exc[
                                "amount"
                            ] * self.carbon_intensity_markets.get(
                                (exc["name"], exc["location"]), {}
                            ).get(
                                "non-fossil", 0.0
                            )

                    if "log parameters" not in ds:
                        ds["log parameters"] = {}

                    # replace current CO2 emissions with new ones
                    if fossil_co2 > 0:
                        for exc in ws.biosphere(
                            ds,
                            ws.equals("name", "Carbon dioxide, fossil"),
                        ):
                            ds["log parameters"]["initial amount of fossil CO2"] = exc[
                                "amount"
                            ]
                            ds["log parameters"]["new amount of fossil CO2"] = float(
                                fossil_co2
                            )
                            exc["amount"] = float(fossil_co2)
                            fossil_co2 = 0

                    if non_fossil_co2 > 0:
                        bio_co2_flows = ws.biosphere(
                            ds,
                            ws.equals("name", "Carbon dioxide, non-fossil"),
                        )

                        for exc in bio_co2_flows:
                            ds["log parameters"]["initial amount of biogenic CO2"] = (
                                exc["amount"]
                            )
                            ds["log parameters"]["new amount of biogenic CO2"] = float(
                                non_fossil_co2
                            )
                            exc["amount"] = float(non_fossil_co2)
                            non_fossil_co2 = 0

                        if non_fossil_co2 > 0 and fossil_co2 == 0:
                            ds["log parameters"]["initial amount of biogenic CO2"] = 0.0
                            ds["log parameters"][
                                "new amount of biogenic CO2"
                            ] = non_fossil_co2

                            ds["exchanges"].append(
                                {
                                    "uncertainty type": 0,
                                    "loc": non_fossil_co2,
                                    "amount": non_fossil_co2,
                                    "name": "Carbon dioxide, non-fossil",
                                    "categories": ("air",),
                                    "type": "biosphere",
                                    "unit": "kilogram",
                                    "tag": "air",
                                }
                            )

                for new_dataset in list(new_ds.values()):
                    self.write_log(new_dataset)
                    # add it to list of created datasets
                    self.add_to_index(new_dataset)
                    self.database.append(new_dataset)

    def write_log(self, dataset, status="created"):
        """
        Write log file.
        """

        logger.info(
            f"{status}|{self.model}|{self.scenario}|{self.year}|"
            f"{dataset['name']}|{dataset['location']}|"
            f"{dataset.get('log parameters', {}).get('initial amount of fossil CO2', '')}|"
            f"{dataset.get('log parameters', {}).get('new amount of fossil CO2', '')}|"
            f"{dataset.get('log parameters', {}).get('initial amount of biogenic CO2', '')}|"
            f"{dataset.get('log parameters', {}).get('new amount of biogenic CO2', '')}"
        )
