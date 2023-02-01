"""
Integrates projections regarding steel production.
"""
import os
from typing import List

from .data_collection import IAMDataCollection
from .transformation import BaseTransformation, ws, wurst
from .utils import DATA_DIR


class Steel(BaseTransformation):
    """
    Class that modifies steel markets in ecoinvent based on IAM output data.

    :ivar database: database dictionary from :attr:`.NewDatabase.database`
    :ivar model: can be 'remind' or 'image'. str from :attr:`.NewDatabase.model`
    :ivar iam_data: xarray that contains IAM data, from :attr:`.NewDatabase.rdc`
    :ivar year: year, from :attr:`.NewDatabase.year`

    """

    def __init__(
        self,
        database: List[dict],
        iam_data: IAMDataCollection,
        model: str,
        pathway: str,
        year: int,
        version: str,
    ) -> None:
        super().__init__(database, iam_data, model, pathway, year)
        self.version = version

    def generate_activities(self):
        """
        This function generates new activities for primary and secondary steel
        production and add them to the wurst database.

        :return: Returns a modified database with newly added steel activities for the corresponding year
        """

        print(f"Log of deleted steel datasets saved in {DATA_DIR / 'logs'}")
        print(f"Log of created steel datasets saved in {DATA_DIR / 'logs'}")

        if not os.path.exists(DATA_DIR / "logs"):
            os.makedirs(DATA_DIR / "logs")

        created_datasets = []

        print("Create steel markets for different regions")

        for i in (
            ("market for steel, low-alloyed", "steel, low-alloyed"),
            ("market for steel, unalloyed", "steel, unalloyed"),
            ("market for steel, chromium steel 18/8", "steel, chromium steel 18/8"),
        ):
            steel_markets = self.fetch_proxies(
                name=i[0],
                ref_prod=i[1],
                production_variable=["steel - primary", "steel - secondary"],
            )

            # adjust share of primary and secondary steel
            if i[0] == "market for steel, low-alloyed":
                for loc, dataset in steel_markets.items():
                    if loc != "World":
                        primary_share = self.iam_data.production_volumes.sel(
                            region=loc, variables="steel - primary"
                        ).interp(year=self.year).values.item(
                            0
                        ) / self.iam_data.production_volumes.sel(
                            region=loc,
                            variables=["steel - primary", "steel - secondary"],
                        ).interp(
                            year=self.year
                        ).sum(
                            dim="variables"
                        ).values.item(
                            0
                        )

                        secondary_share = 1 - primary_share

                        new_exc = [
                            {
                                "uncertainty type": 0,
                                "loc": primary_share,
                                "amount": primary_share,
                                "type": "technosphere",
                                "production volume": 1,
                                "product": "steel, low-alloyed",
                                "name": "steel production, converter, low-alloyed",
                                "unit": "kilogram",
                                "location": loc,
                            },
                            {
                                "uncertainty type": 0,
                                "loc": secondary_share,
                                "amount": secondary_share,
                                "type": "technosphere",
                                "production volume": 1,
                                "product": "steel, low-alloyed",
                                "name": "steel production, electric, low-alloyed",
                                "unit": "kilogram",
                                "location": loc,
                            },
                        ]

                        dataset["exchanges"] = [
                            e
                            for e in dataset["exchanges"]
                            if e["type"] == "production" or e["unit"] == "ton kilometer"
                        ]
                        dataset["exchanges"].extend(new_exc)

            else:
                for loc, dataset in steel_markets.items():
                    if loc != "World":
                        name_ref = [
                            (e["name"], e.get("product"))
                            for e in dataset["exchanges"]
                            if "steel production" in e["name"]
                        ][0]
                        name, ref = name_ref

                        dataset["exchanges"] = [
                            e
                            for e in dataset["exchanges"]
                            if e["type"] == "production" or e["unit"] == "ton kilometer"
                        ]

                        dataset["exchanges"].append(
                            {
                                "uncertainty type": 0,
                                "loc": 1,
                                "amount": 1,
                                "type": "technosphere",
                                "production volume": 1,
                                "product": ref,
                                "name": name,
                                "unit": "kilogram",
                                "location": loc,
                            }
                        )

            # populate World dataset
            steel_markets["World"]["exchanges"] = [
                x
                for x in steel_markets["World"]["exchanges"]
                if x["type"] == "production"
            ]
            regions = [r for r in self.regions if r != "World"]

            for region in regions:
                share = (
                    self.iam_data.production_volumes.sel(
                        variables=["steel - primary", "steel - secondary"],
                        region=region,
                    )
                    .interp(year=self.year)
                    .sum(dim="variables")
                    / self.iam_data.production_volumes.sel(
                        variables=["steel - primary", "steel - secondary"],
                        region="World",
                    )
                    .interp(year=self.year)
                    .sum(dim="variables")
                ).values.item(0)

                steel_markets["World"]["exchanges"].append(
                    {
                        "name": i[0],
                        "product": i[1],
                        "amount": share,
                        "unit": "kilogram",
                        "type": "technosphere",
                        "location": region,
                    }
                )

            self.database.extend(list(steel_markets.values()))

            created_datasets.extend(
                [
                    (act["name"], act["reference product"], act["location"])
                    for act in steel_markets.values()
                ]
            )

        # Determine all steel activities in the database. Empty old datasets.
        print("Create new steel production datasets and empty old datasets")
        d_act_primary_steel = {
            mat: self.fetch_proxies(
                name=mat[0],
                ref_prod=mat[1],
                production_variable=["steel - primary"],
                relink=True,
            )
            for mat in zip(
                self.material_map["steel - primary"],
                ["steel"] * len(self.material_map["steel - primary"]),
            )
        }
        d_act_secondary_steel = {
            mat: self.fetch_proxies(
                name=mat[0],
                ref_prod=mat[1],
                production_variable=["steel - secondary"],
                relink=True,
            )
            for mat in zip(
                self.material_map["steel - secondary"],
                ["steel"] * len(self.material_map["steel - secondary"]),
            )
        }

        # Create region-specific pig iron production datasets

        print("Create pig iron production datasets")

        pig_iron_production = {
            mat: self.fetch_proxies(
                name=mat[0], ref_prod=mat[1], production_variable=["steel - primary"]
            )
            for mat in zip(
                self.material_map["pig iron"],
                ["pig iron"] * len(self.material_map["pig iron"]),
            )
        }

        for _, dataset in pig_iron_production.items():
            self.database.extend(list(dataset.values()))

        pig_iron_markets = self.fetch_proxies(
            name="market for pig iron",
            ref_prod="pig iron",
            production_variable=["steel - primary"],
        )
        self.database.extend(list(pig_iron_markets.values()))

        d_act_steel = {
            **pig_iron_production,
            **d_act_primary_steel,
            **d_act_secondary_steel,
        }

        # Scale down fuel exchanges, according to efficiency improvement as
        # forecast by the IAM:
        list_fuels = [
            "diesel",
            "coal",
            "lignite",
            "coke",
            "fuel",
            "meat",
            "gas",
            "oil",
            "electricity",
            "natural gas",
            "steam",
        ]

        for steel in d_act_steel:
            for region, activity in d_act_steel[steel].items():
                # the correction factor applied to all fuel/electricity input is
                # equal to the ration fuel/output in the year in question
                # divided by the ratio fuel/output in 2020

                sector = (
                    "steel - primary"
                    if any(i in activity["name"] for i in ["converter", "pig iron"])
                    else "steel - secondary"
                )
                scaling_factor = 1 / self.find_iam_efficiency_change(
                    variable=sector, location=activity["location"]
                )

                # update comments
                text = (
                    "This dataset has been modified by `premise`, according to "
                    f"the performance for steel production indicated by the IAM model {self.model.upper()} "
                    f"for the IAM region {region} in {self.year}, following the scenario {self.scenario}. "
                    f"The energy efficiency of the process "
                    f"has been improved by {int((1 - scaling_factor) * 100)}%. "
                )

                d_act_steel[steel][region]["comment"] = text + activity["comment"]

                wurst.change_exchanges_by_constant_factor(
                    activity,
                    scaling_factor,
                    technosphere_filters=[
                        ws.either(*[ws.contains("name", x) for x in list_fuels])
                    ],
                    biosphere_filters=[ws.contains("name", "Carbon dioxide, fossil")],
                )

                # Add carbon capture-related energy exchanges
                # Carbon capture rate: share of capture of total CO2 emitted
                # Note: only if variables exist in IAM data

                # Carbon capture rate: share of capture of total CO2 emitted
                carbon_capture_rate = self.get_carbon_capture_rate(
                    loc=activity["location"], sector="steel"
                )

                if carbon_capture_rate > 0:
                    for co2_flow in ws.biosphere(
                        activity, ws.contains("name", "Carbon dioxide, fossil")
                    ):
                        co2_amount = co2_flow["amount"]
                        co2_emitted = co2_amount * (1 - carbon_capture_rate)
                        co2_flow["amount"] = co2_emitted

                        # create the CCS dataset to fit this clinker production dataset
                        # and add it to the database
                        self.create_ccs_dataset(
                            region,
                            bio_co2_stored=0,
                            bio_co2_leaked=0,
                        )

                        # add an input from this CCS dataset in the clinker dataset
                        ccs_exc = {
                            "uncertainty type": 0,
                            "loc": 0,
                            "amount": co2_amount - co2_emitted,
                            "type": "technosphere",
                            "production volume": 0,
                            "name": "CO2 capture, at cement production plant, "
                            "with underground storage, post, 200 km",
                            "unit": "kilogram",
                            "location": activity["location"],
                            "product": "CO2, captured and stored",
                        }
                        activity["exchanges"].append(ccs_exc)

                # Update hot pollutant emission according to GAINS
                self.update_pollutant_emissions(dataset=activity, sector="steel")

            if steel != ("pig iron production", "pig iron"):
                self.database.extend(list(d_act_steel[steel].values()))

        print("Done!")
