"""
cement.py contains the class `Cement`, which inherits from `BaseTransformation`.
This class transforms the cement markets and clinker and cement production activities
of the wurst database, based on projections from the IAM scenario.
It also the generic market for cement to reflect the projected clinker-to-cement ratio.
It eventually re-links all the cement-consuming activities (e.g., concrete production)
of the wurst database to the newly created cement markets.

"""
import csv
import os
from datetime import date
from pathlib import Path

from .transformation import (
    BaseTransformation,
    Dict,
    IAMDataCollection,
    List,
    get_shares_from_production_volume,
    get_suppliers_of_a_region,
    np,
    remove_exchanges,
    ws,
)
from .utils import DATA_DIR, get_clinker_ratio_ecoinvent, get_clinker_ratio_remind


class Cement(BaseTransformation):
    """
    Class that modifies clinker and cement production datasets in ecoinvent.
    It creates region-specific new clinker production datasets (and deletes the original ones).
    It adjusts the kiln efficiency based on the improvement indicated in the IAM file, relative to 2020.
    It adjusts non-CO2 pollutants emission, based on improvement indicated by the GAINS file, relative to 2020.
    It adds CCS, if indicated in the IAM file.
    It creates regions-specific cement production datasets (and deletes the original ones).
    It adjust electricity consumption in cement production datasets.
    It creates regions-specific cement market datasets (and deletes the original ones).
    It adjust the clinker-to-cement ratio in the generic cement market dataset.


    :ivar database: wurst database, which is a list of dictionaries
    :ivar iam_data: IAM data
    :ivar model: name of the IAM model (e.g., "remind", "image")
    :ivar pathway: name of the IAM pathway (e.g., "SSP2-Base")
    :ivar year: year of the pathway (e.g., 2030)
    :ivar version: version of ecoinvent database (e.g., "3.7")

    """

    def __init__(
        self,
        database: List[dict],
        iam_data: IAMDataCollection,
        model: str,
        pathway: str,
        year: int,
        version: str,
    ):
        super().__init__(database, iam_data, model, pathway, year)
        self.version = version
        self.clinker_ratio_eco = get_clinker_ratio_ecoinvent(version)
        self.clinker_ratio_remind = get_clinker_ratio_remind(self.year)

    def build_clinker_production_datasets(self) -> Dict[str, dict]:
        """
        Builds clinker production datasets for each IAM region.
        Adds CO2 capture and Storage, if needed.

        :return: a dictionary with IAM regions as keys and clinker production datasets as values.
        """

        # Fetch clinker production activities and store them in a dictionary
        d_act_clinker = self.fetch_proxies(
            name="clinker production",
            ref_prod="clinker",
            production_variable="cement",
            relink=True,
        )

        # Fuel exchanges to remove
        list_fuels = [
            "diesel",
            "coal",
            "lignite",
            "coke",
            "fuel",
            "meat",
            "gas",
            "oil",
            "wood",
            "waste",
        ]

        # we first create current clinker production activities for each region
        # from GNR data, to get a new fuel efficiency and mix than what is
        # currently in ecoinvent

        # Remove fuel and electricity exchanges in each activity
        d_act_clinker = remove_exchanges(d_act_clinker, list_fuels)

        for region, dataset in d_act_clinker.items():
            # Production volume by kiln type
            energy_input_per_kiln_type = self.iam_data.gnr_data.sel(
                region=self.geo.iam_to_iam_region(region, from_iam="remind")
                if self.model == "image"
                else region,
                variables=[
                    v
                    for v in self.iam_data.gnr_data.variables.values
                    if "Production volume share" in v
                ],
            ).clip(0, 1)

            # Energy input per ton of clinker, in MJ, per kiln type
            energy_input_per_kiln_type /= energy_input_per_kiln_type.sum(axis=0)

            energy_eff_per_kiln_type = self.iam_data.gnr_data.sel(
                region=self.geo.iam_to_iam_region(region, from_iam="remind")
                if self.model == "image"
                else region,
                variables=[
                    v
                    for v in self.iam_data.gnr_data.variables.values
                    if "Thermal energy consumption" in v
                ],
            )

            # Weighted average energy input per ton clinker, in MJ
            energy_input_per_ton_clinker = (
                energy_input_per_kiln_type.values * energy_eff_per_kiln_type.values
            )

            # the correction factor applied to all fuel/electricity input is
            # equal to the ratio fuel/output in the year in question
            # divided by the ratio fuel/output in 2020

            scaling_factor = 1 / self.find_iam_efficiency_change(
                variable="cement", location=dataset["location"]
            )
            energy_input_per_ton_clinker = energy_input_per_ton_clinker.sum()
            energy_input_per_ton_clinker *= scaling_factor

            # Limit the energy input to the upper bound value
            # of the theoretical minimum considered by the IEA of
            # 2.8 GJ/t clinker.
            energy_input_per_ton_clinker = np.clip(
                energy_input_per_ton_clinker, 2800, None
            )

            # Fuel mix (waste, biomass, fossil)
            fuel_mix = self.iam_data.gnr_data.sel(
                variables=[
                    "Share waste fuel",
                    "Share biomass fuel",
                    "Share fossil fuel",
                ],
                region=self.geo.iam_to_iam_region(region, from_iam="remind")
                if self.model == "image"
                else region,
            ).clip(0, 1)

            fuel_mix /= fuel_mix.sum(axis=0)

            # Calculate quantities (in kg) of fuel, per type of fuel, per ton of clinker
            # MJ per ton of clinker * fuel mix * (1 / lower heating value)
            fuel_qty_per_type = (
                energy_input_per_ton_clinker
                * fuel_mix
                * 1
                / np.array(
                    [
                        float(self.fuels_specs["waste"]["lhv"]),
                        float(self.fuels_specs["wood pellet"]["lhv"]),
                        float(self.fuels_specs["hard coal"]["lhv"]),
                    ]
                )
            )

            fuel_fossil_co2_per_type = (
                energy_input_per_ton_clinker
                * fuel_mix
                * np.array(
                    [
                        (
                            self.fuels_specs["waste"]["co2"]
                            * (1 - self.fuels_specs["waste"]["biogenic_share"])
                        ),
                        (
                            self.fuels_specs["wood pellet"]["co2"]
                            * (1 - self.fuels_specs["wood pellet"]["biogenic_share"])
                        ),
                        (
                            self.fuels_specs["hard coal"]["co2"]
                            * (1 - self.fuels_specs["hard coal"]["biogenic_share"])
                        ),
                    ]
                )
            )

            fuel_biogenic_co2_per_type = (
                energy_input_per_ton_clinker.sum()
                * fuel_mix
                * np.array(
                    [
                        (
                            self.fuels_specs["waste"]["co2"]
                            * (self.fuels_specs["waste"]["biogenic_share"])
                        ),
                        (
                            self.fuels_specs["wood pellet"]["co2"]
                            * (self.fuels_specs["wood pellet"]["biogenic_share"])
                        ),
                        (
                            self.fuels_specs["hard coal"]["co2"]
                            * (self.fuels_specs["hard coal"]["biogenic_share"])
                        ),
                    ]
                )
            )

            # Append it to the dataset exchanges
            new_exchanges = []

            for f_idx, fuel in enumerate(
                [
                    ("waste", "waste plastic, mixture"),
                    ("wood pellet", "wood pellet, measured as dry mass"),
                    ("hard coal", "hard coal"),
                ]
            ):
                # Select waste fuel providers, fitting the IAM region
                # Fetch respective shares based on production volumes
                # if they cannot be found for the ecoinvent locations concerned
                # we widen the scope to EU-based datasets, and RoW
                ecoinvent_regions = self.geo.iam_to_ecoinvent_location(region)
                possible_locations = [
                    ecoinvent_regions,
                    ["RER", "Europe without Switzerland"],
                    ["RoW", "GLO"],
                ]
                suppliers, counter = [], 0

                while len(suppliers) == 0:
                    suppliers = list(
                        get_suppliers_of_a_region(
                            database=self.database,
                            locations=possible_locations[counter],
                            names=list(self.fuel_map[fuel[0]]),
                            reference_product=fuel[1],
                            unit="kilogram",
                            exclude=["ash", "mine"],
                        )
                    )
                    counter += 1

                suppliers = get_shares_from_production_volume(suppliers)

                for supplier, share in suppliers.items():
                    new_exchanges.append(
                        {
                            "uncertainty type": 0,
                            "loc": 1,
                            "amount": (share * fuel_qty_per_type[f_idx].values) / 1000,
                            "type": "technosphere",
                            "production volume": 0,
                            "product": supplier[2],
                            "name": supplier[0],
                            "unit": supplier[-1],
                            "location": supplier[1],
                        }
                    )

            dataset["exchanges"].extend(new_exchanges)

            dataset["exchanges"] = [v for v in dataset["exchanges"] if v]

            # Carbon capture rate: share of capture of total CO2 emitted
            carbon_capture_rate = self.get_carbon_capture_rate(
                loc=dataset["location"], sector="cement"
            )

            # Update fossil CO2 exchange,
            # add 525 kg of fossil CO_2 from calcination
            co2_exc = [
                e for e in dataset["exchanges"] if e["name"] == "Carbon dioxide, fossil"
            ]

            if co2_exc:
                fossil_co2_exc = co2_exc[0]
                fossil_co2_exc["amount"] = (
                    (fuel_fossil_co2_per_type.sum().values + 525) / 1000
                ) * (1 - carbon_capture_rate)
                fossil_co2_exc["uncertainty type"] = 0

            else:
                # the fossil CO2 flow does not exist
                amount = ((fuel_fossil_co2_per_type.sum().values + 525) / 1000) * (
                    1 - carbon_capture_rate
                )
                fossil_co2_exc = {
                    "uncertainty type": 0,
                    "loc": float(amount),
                    "amount": float(amount),
                    "type": "biosphere",
                    "name": "Carbon dioxide, fossil",
                    "unit": "kilogram",
                    "categories": ("air",),
                }
                dataset["exchanges"].append(fossil_co2_exc)

            co2_exc = [
                e
                for e in dataset["exchanges"]
                if e["name"] == "Carbon dioxide, non-fossil"
            ]

            if co2_exc:
                # Update biogenic CO2 exchange
                biogenic_co2_exc = co2_exc[0]
                biogenic_co2_exc["amount"] = (
                    fuel_biogenic_co2_per_type.sum().values / 1000
                ) * (1 - carbon_capture_rate)
                biogenic_co2_exc["uncertainty type"] = 0

            else:
                # There isn't a biogenic CO2 emission exchange
                amount = (fuel_biogenic_co2_per_type.sum().values / 1000) * (
                    1 - carbon_capture_rate
                )
                biogenic_co2_exc = {
                    "uncertainty type": 0,
                    "loc": float(amount),
                    "amount": float(amount),
                    "type": "biosphere",
                    "name": "Carbon dioxide, non-fossil",
                    "unit": "kilogram",
                    "input": ("biosphere3", "eba59fd6-f37e-41dc-9ca3-c7ea22d602c7"),
                    "categories": ("air",),
                }
                dataset["exchanges"].append(biogenic_co2_exc)

            # add CCS-related dataset
            if carbon_capture_rate > 0:
                # total CO2 emissions = bio CO2 emissions
                # + fossil CO2 emissions
                # + calcination emissions
                total_co2_emissions = (
                    fuel_fossil_co2_per_type.sum()
                    + fuel_biogenic_co2_per_type.sum()
                    + 525
                )
                # share bio CO2 stored = sum of biogenic fuel emissions / total CO2 emissions
                bio_co2_stored = (
                    fuel_biogenic_co2_per_type.sum() / total_co2_emissions
                ).values.item(0)

                # 0.11 kg CO2 leaks per kg captured
                # we need to align the CO2 composition with
                # the CO2 composition of the cement plant
                bio_co2_leaked = (
                    fuel_biogenic_co2_per_type.sum() / total_co2_emissions
                ).values.item(0) * 0.11

                # create the CCS dataset to fit this clinker production dataset
                # and add it to the database
                self.create_ccs_dataset(
                    region,
                    bio_co2_stored,
                    bio_co2_leaked,
                )

                # add an input from this CCS dataset in the clinker dataset
                ccs_exc = {
                    "uncertainty type": 0,
                    "loc": 0,
                    "amount": float((total_co2_emissions / 1000) * carbon_capture_rate),
                    "type": "technosphere",
                    "production volume": 0,
                    "name": "CO2 capture, at cement production plant, with underground storage, post, 200 km",
                    "unit": "kilogram",
                    "location": dataset["location"],
                    "product": "CO2, captured and stored",
                }
                dataset["exchanges"].append(ccs_exc)

            dataset["exchanges"] = [v for v in dataset["exchanges"] if v]

            share_fossil_from_fuel = int(
                (
                    fuel_fossil_co2_per_type.sum()
                    / (fuel_fossil_co2_per_type.sum() + 525)
                )
                * 100
            )

            share_fossil_from_calcination = 100 - int(
                (
                    fuel_fossil_co2_per_type.sum()
                    / np.sum(fuel_fossil_co2_per_type.sum() + 525)
                )
                * 100
            )

            dataset["comment"] = (
                "Dataset modified by `premise` based on WBCSD's GNR data and IAM projections "
                + " for the cement industry.\n"
                + f"Calculated energy input per kg clinker: {np.round(energy_input_per_ton_clinker.sum(), 1) / 1000}"
                f" MJ/kg clinker.\n"
                + f"Improvement of specific energy use compared to 2020: {(scaling_factor - 1) * 100} %.\n"
                + f"Share of biomass fuel energy-wise: {int(fuel_mix[1] * 100)} pct.\n"
                + f"Share of waste fuel energy-wise: {int(fuel_mix[0] * 100)} pct.\n"
                + f"Share of biogenic carbon in waste fuel energy-wise: {int(self.fuels_specs['waste']['biogenic_share'] * 100)} pct.\n"
                + f"Share of fossil CO2 emissions from fuel combustion: {share_fossil_from_fuel} pct.\n"
                + f"Share of fossil CO2 emissions from calcination: {share_fossil_from_calcination} pct.\n"
                + f"Rate of carbon capture: {int(carbon_capture_rate * 100)} pct.\n"
            ) + dataset["comment"]

        print(
            "Adjusting emissions of hot pollutants for clinker production datasets..."
        )
        d_act_clinker = {
            k: self.update_pollutant_emissions(dataset=v, sector="cement")
            for k, v in d_act_clinker.items()
        }

        return d_act_clinker

    def adjust_clinker_ratio(self, d_act: Dict[str, dict]) -> Dict[str, dict]:
        """Adjust the cement suppliers composition for "cement, unspecified", in order to reach
        the average clinker-to-cement ratio given by the IAM.

        The supply of the cement with the highest clinker-to-cement ratio is decreased by 1% to the favor of
        the supply of the cement with the lowest clinker-to-cement ratio, and the average clinker-to-cement ratio
        is calculated.

        This operation is repeated until the average clinker-to-cement ratio aligns with that given by the IAM.
        When the supply of the cement with the highest clinker-to-cement ratio goes below 1%,
        the cement with the second highest clinker-to-cement ratio becomes affected and so forth.

        """

        for region in d_act:
            # we fetch the clinker-to-cement ratio forecast
            # by REMIND (as other IAMs do not seem to consider
            # this). This ratio decreases over time
            # (reducing the amount of clinker in cement), but at
            # different pace across regions.
            ratio_to_reach = self.clinker_ratio_remind.sel(
                dict(
                    region=self.geo.iam_to_iam_region(region, from_iam="remind")
                    if self.model == "image"
                    else region
                )
            ).values

            share = []
            ratio = []

            # we loop through the exchange of the dataset
            # if we meet a cement exchange, we store its clinker-to-cement ratio
            # in `ratio`
            for exc in d_act[region]["exchanges"]:
                if "cement" in exc["product"] and exc["type"] == "technosphere":
                    share.append(exc["amount"])
                    try:
                        ratio.append(
                            self.clinker_ratio_eco[(exc["name"], exc["location"])]
                        )
                    except KeyError:
                        print(
                            f"clinker-to-cement ratio not found for {exc['name'], exc['location']}"
                        )

            share = np.array(share)
            ratio = np.array(ratio)

            # once we know what cement exchanges are in the dataset
            # we can calculate the volume-weighted clinker-to-cement ratio
            # of the dataset
            average_ratio = (share * ratio).sum()

            # As long as we do not reach the clinker-to-cement ratio forecast by REMIND
            # we will decrease the input of the highest clinker-containing cement
            # and increase the input of the lowest clinker-containing cement
            iteration = 0
            while average_ratio > ratio_to_reach and iteration < 100:
                share[share == 0] = np.nan

                ratio = np.where(share >= 0.001, ratio, np.nan)

                highest_ratio = np.nanargmax(ratio)
                lowest_ratio = np.nanargmin(ratio)

                share[highest_ratio] -= 0.01
                share[lowest_ratio] += 0.01

                average_ratio = (np.nan_to_num(ratio) * np.nan_to_num(share)).sum()
                iteration += 1

            share = np.nan_to_num(share)

            count = 0
            for exc in d_act[region]["exchanges"]:
                if "cement" in exc["product"] and exc["type"] == "technosphere":
                    exc["amount"] = share[count]
                    count += 1

        return d_act

    def update_electricity_exchanges(self, d_act: Dict[str, dict]) -> Dict[str, dict]:
        """
        Update electricity exchanges in cement production datasets.
        Here, we use data from the GNR database for current consumption, and the 2018 IEA/GNR roadmap publication
        for future consumption.
        Electricity consumption equals electricity use minus on-site electricity generation from excess heat recovery.

        :return:
        """
        d_act = remove_exchanges(d_act, ["electricity"])

        for region in d_act:
            new_exchanges = []

            electricity = (
                self.iam_data.gnr_data.loc[
                    dict(
                        variables=["Power generation", "Power consumption"],
                        region=self.geo.iam_to_gains_region(region),
                    )
                ]
                / 1000
            )
            electricity_needed = (
                electricity.sel(variables="Power consumption")
                - electricity.sel(variables="Power generation")
            ).values

            # Fetch electricity-producing technologies contained in the IAM region
            # if they cannot be found for the ecoinvent locations concerned
            # we widen the scope to EU-based datasets, and RoW
            ecoinvent_regions = self.geo.iam_to_ecoinvent_location(region)
            possible_locations = [[region], ecoinvent_regions, ["RER"], ["RoW"]]
            suppliers, counter = [], 0

            while len(suppliers) == 0:
                suppliers = list(
                    get_suppliers_of_a_region(
                        database=self.database,
                        locations=possible_locations[counter],
                        names=["electricity, medium voltage"],
                        reference_product="electricity",
                        unit="kilowatt hour",
                    )
                )
                counter += 1

            suppliers = get_shares_from_production_volume(suppliers)

            for supplier, share in suppliers.items():
                new_exchanges.append(
                    {
                        "uncertainty type": 0,
                        "loc": electricity_needed * share,
                        "amount": electricity_needed * share,
                        "type": "technosphere",
                        "production volume": 0,
                        "product": supplier[2],
                        "name": supplier[0],
                        "unit": supplier[-1],
                        "location": supplier[1],
                    }
                )

            d_act[region]["exchanges"].extend(new_exchanges)
            d_act[region]["exchanges"] = [v for v in d_act[region]["exchanges"] if v]

            txt = (
                "Dataset modified by `premise` based on WBCSD's GNR data and 2018 IEA "
                "roadmap for the cement industry.\n "
                f"Electricity consumption per kg cement: {electricity_needed} kWh."
                f"Of which a part was generated from on-site waste heat recovery.\n"
            )

            if "comment" in d_act[region]:
                d_act[region]["comment"] += txt
            else:
                d_act[region]["comment"] = txt

        return d_act

    def add_datasets_to_database(self) -> None:
        """
        Runs a series of methods that create new clinker and cement production datasets
        and new cement market datasets.
        :return: Does not return anything. Modifies in place.
        """

        print("\nStart integration of cement data...\n")

        print(f"Log of deleted cement datasets saved in {DATA_DIR}/logs")
        print(f"Log of created cement datasets saved in {DATA_DIR}/logs")

        if not os.path.exists(DATA_DIR / "logs"):
            os.makedirs(DATA_DIR / "logs")

        with open(
            DATA_DIR
            / f"logs/log deleted cement datasets {self.model} {self.scenario} {self.year}-{date.today()}.csv",
            "w",
            encoding="utf-8",
        ) as csv_file:
            writer = csv.writer(csv_file, delimiter=";", lineterminator="\n")
            writer.writerow(["dataset name", "reference product", "location"])

        with open(
            DATA_DIR
            / f"logs/log created cement datasets {self.model} {self.scenario} {self.year}-{date.today()}.csv",
            "w",
            encoding="utf-8",
        ) as csv_file:
            writer = csv.writer(csv_file, delimiter=";", lineterminator="\n")
            writer.writerow(["dataset name", "reference product", "location"])

        created_datasets = []

        print("\nCreate new clinker production datasets and delete old datasets")
        clinker_prod_datasets = list(self.build_clinker_production_datasets().values())
        self.database.extend(clinker_prod_datasets)

        created_datasets.extend(
            [
                (act["name"], act["reference product"], act["location"])
                for act in clinker_prod_datasets
            ]
        )

        print("\nCreate new clinker market datasets and delete old datasets")
        clinker_market_datasets = list(
            self.fetch_proxies(
                name="market for clinker",
                ref_prod="clinker",
                production_variable="cement",
            ).values()
        )

        self.database.extend(clinker_market_datasets)

        created_datasets.extend(
            [
                (act["name"], act["reference product"], act["location"])
                for act in clinker_market_datasets
            ]
        )

        print('Adjust clinker-to-cement ratio in "unspecified cement" datasets')

        if self.version == 3.5:
            name = "market for cement, unspecified"
            ref_prod = "cement, unspecified"

        else:
            name = "cement, all types to generic market for cement, unspecified"
            ref_prod = "cement, unspecified"

        act_cement_unspecified = self.fetch_proxies(
            name=name,
            ref_prod=ref_prod,
            production_variable="cement",
        )
        act_cement_unspecified = self.adjust_clinker_ratio(act_cement_unspecified)
        self.database.extend(list(act_cement_unspecified.values()))

        created_datasets.extend(
            [
                (act["name"], act["reference product"], act["location"])
                for act in act_cement_unspecified.values()
            ]
        )

        print("\nCreate new cement market datasets")

        # cement markets
        markets = ws.get_many(
            self.database,
            ws.contains("name", "market for cement"),
            ws.contains("reference product", "cement"),
            ws.doesnt_contain_any("name", ["factory", "tile", "sulphate", "plaster"]),
            ws.doesnt_contain_any("location", self.regions),
        )

        markets = {(m["name"], m["reference product"]) for m in markets}

        for dataset in markets:
            new_cement_markets = self.fetch_proxies(
                name=dataset[0],
                ref_prod=dataset[1],
                production_variable="cement",
            )

            self.database.extend(list(new_cement_markets.values()))
            created_datasets.extend(
                [
                    (act["name"], act["reference product"], act["location"])
                    for act in new_cement_markets.values()
                ]
            )

        print(
            "\nCreate new cement production datasets and adjust electricity consumption"
        )
        # cement production
        production = ws.get_many(
            self.database,
            ws.contains("name", "cement production"),
            ws.contains("reference product", "cement"),
            ws.doesnt_contain_any("name", ["factory", "tile", "sulphate", "plaster"]),
        )

        production = {(p["name"], p["reference product"]) for p in production}

        for dataset in production:
            # Fetch proxy datasets (one per IAM region)
            # Delete old datasets
            new_cement_production = self.fetch_proxies(
                name=dataset[0],
                ref_prod=dataset[1],
                production_variable="cement",
            )
            # Update electricity use
            new_cement_production = self.update_electricity_exchanges(
                new_cement_production
            )

            # add them to the wurst database
            self.database.extend(list(new_cement_production.values()))

            # relink cement production-consuming datasets to the new ones
            created_datasets.extend(
                [
                    (act["name"], act["reference product"], act["location"])
                    for act in new_cement_production.values()
                ]
            )

        Path(DATA_DIR / "logs").mkdir(parents=True, exist_ok=True)

        with open(
            DATA_DIR
            / f"logs/log created cement datasets {self.model} {self.scenario} {self.year}-{date.today()}.csv",
            "a",
            encoding="utf-8",
        ) as csv_file:
            writer = csv.writer(csv_file, delimiter=";", lineterminator="\n")
            for line in created_datasets:
                writer.writerow(line)

        print("Done!")
