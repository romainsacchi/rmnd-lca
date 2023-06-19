"""
cement.py contains the class `Cement`, which inherits from `BaseTransformation`.
This class transforms the cement markets and clinker and cement production activities
of the wurst database, based on projections from the IAM scenario.
It also the generic market for cement to reflect the projected clinker-to-cement ratio.
It eventually re-links all the cement-consuming activities (e.g., concrete production)
of the wurst database to the newly created cement markets.

"""

import logging.config
from collections import defaultdict
from pathlib import Path

import yaml

from .transformation import BaseTransformation, Dict, IAMDataCollection, List, np, ws
from .utils import DATA_DIR

LOG_CONFIG = DATA_DIR / "utils" / "logging" / "logconfig.yaml"
# directory for log files
DIR_LOG_REPORT = Path.cwd() / "export" / "logs"
# if DIR_LOG_REPORT folder does not exist
# we create it
if not Path(DIR_LOG_REPORT).exists():
    Path(DIR_LOG_REPORT).mkdir(parents=True, exist_ok=True)

with open(LOG_CONFIG, "r") as f:
    config = yaml.safe_load(f.read())
    logging.config.dictConfig(config)

logger = logging.getLogger("cement")


class Cement(BaseTransformation):
    """
    Class that modifies clinker and cement production datasets in ecoinvent.
    It creates region-specific new clinker production datasets (and deletes the original ones).
    It adjusts the kiln efficiency based on the improvement indicated in the IAM file, relative to 2020.
    It adds CCS, if indicated in the IAM file.
    It creates regions-specific cement production datasets (and deletes the original ones).
    It adjusts electricity consumption in cement production datasets.
    It creates regions-specific cement market datasets (and deletes the original ones).
    It adjusts the clinker-to-cement ratio in the generic cement market dataset.


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
        system_model: str,
        modified_datasets: dict,
    ):
        super().__init__(
            database,
            iam_data,
            model,
            pathway,
            year,
            version,
            system_model,
            modified_datasets,
        )
        self.version = version

    def fetch_current_energy_details(self, dataset):
        """
        Fetches the current energy consumption per ton of clinker

        :param dataset: clinker production dataset
        :return: current energy consumption per ton of clinker, in kj/kg clinker
        """

        # store in a dictionary the amount of fuel and electricity input
        # along with the energy and fossil and biogenic CO2
        # emissions associated with the combustion of each fuel

        _ = lambda x: self.calculate_input_energy(x["name"], x["amount"], x["unit"])

        d_fuels = defaultdict(dict)

        for exc in dataset["exchanges"]:
            if (
                exc["name"] in self.cement_fuels_map["cement, dry feed rotary kiln"]
                and exc["type"] == "technosphere"
            ):
                if exc["name"] not in d_fuels:
                    d_fuels[exc["name"]] = {
                        "amount": exc["amount"],
                        "energy": _(exc) * 1000,
                        "fossil CO2": self.fuels_specs[
                            self.fuel_map_reverse[exc["name"]]
                        ]["co2"]
                        * _(exc)
                        * (
                            1
                            - self.fuels_specs[self.fuel_map_reverse[exc["name"]]][
                                "biogenic_share"
                            ]
                        ),
                        "biogenic CO2": self.fuels_specs[
                            self.fuel_map_reverse[exc["name"]]
                        ]["co2"]
                        * _(exc)
                        * self.fuels_specs[self.fuel_map_reverse[exc["name"]]][
                            "biogenic_share"
                        ],
                    }
                else:
                    d_fuels[exc["name"]]["amount"] += exc["amount"]
                    d_fuels[exc["name"]]["energy"] += _(exc) * 1000
                    d_fuels[exc["name"]]["fossil CO2"] += (
                        self.fuels_specs[self.fuel_map_reverse[exc["name"]]]["co2"]
                        * _(exc)
                        * (
                            1
                            - self.fuels_specs[self.fuel_map_reverse[exc["name"]]][
                                "biogenic_share"
                            ]
                        )
                    )
                    d_fuels[exc["name"]]["biogenic CO2"] += (
                        self.fuels_specs[self.fuel_map_reverse[exc["name"]]]["co2"]
                        * _(exc)
                        * self.fuels_specs[self.fuel_map_reverse[exc["name"]]][
                            "biogenic_share"
                        ]
                    )

        return d_fuels

    def fetch_bio_ratio(self, dataset: dict) -> float:
        """
        Return teh ratio between bio and non-bio CO2 emissions
        in the dataset.

        :param dataset: clinker production dataset
        :return: ratio between bio and non-bio CO2 emissions
        """

        bio_co2 = sum(
            [
                e["amount"]
                for e in dataset["exchanges"]
                if e["name"] == "Carbon dioxide, non-fossil"
            ]
        )
        non_bio_co2 = sum(
            [
                e["amount"]
                for e in dataset["exchanges"]
                if e["name"] == "Carbon dioxide, fossil"
            ]
        )

        return bio_co2 / (bio_co2 + non_bio_co2)

    def rescale_fuel_inputs(self, dataset, scaling_factor, energy_details):
        if scaling_factor != 1:
            for exc in dataset["exchanges"]:
                if exc["name"] in self.cement_fuels_map["cement, dry feed rotary kiln"]:
                    exc["amount"] *= scaling_factor

                    # update energy_details
                    energy_details[exc["name"]]["amount"] *= scaling_factor
                    energy_details[exc["name"]]["energy"] *= scaling_factor
                    energy_details[exc["name"]]["fossil CO2"] *= scaling_factor
                    energy_details[exc["name"]]["biogenic CO2"] *= scaling_factor

        new_energy_input = sum([d["energy"] for d in energy_details.values()])

        dataset["log parameters"].update(
            {
                "new energy input per ton clinker": new_energy_input,
            }
        )

        return dataset

    def rescale_emissions(self, dataset, energy_details, scaling_factor):
        if scaling_factor != 1:
            for exc in dataset["exchanges"]:
                if exc["name"].lower().startswith("carbon dioxide"):
                    if "non-fossil" in exc["name"].lower():
                        dataset["log parameters"].update(
                            {
                                "initial biogenic CO2": exc["amount"],
                            }
                        )
                        # calculate total biogenic CO2 from energy_details
                        total_biogenic_CO2 = sum(
                            [d["biogenic CO2"] for d in energy_details.values()]
                        )

                        exc["amount"] = total_biogenic_CO2
                        dataset["log parameters"].update(
                            {
                                "new biogenic CO2": total_biogenic_CO2,
                            }
                        )
                    else:
                        # calculate total fossil CO2 from energy_details
                        total_fossil_CO2 = sum(
                            [d["fossil CO2"] for d in energy_details.values()]
                        )

                        dataset["log parameters"].update(
                            {
                                "initial fossil CO2": exc["amount"],
                            }
                        )

                        # remove 525 kg for calcination
                        exc["amount"] = 0.525 + total_fossil_CO2
                        dataset["log parameters"].update(
                            {
                                "new fossil CO2": exc["amount"],
                            }
                        )
                else:
                    if exc["type"] == "biosphere":
                        exc["amount"] *= scaling_factor

        else:
            for exc in dataset["exchanges"]:
                if exc["name"].lower().startswith("carbon dioxide"):
                    if "non-fossil" in exc["name"].lower():
                        dataset["log parameters"].update(
                            {
                                "initial biogenic CO2": exc["amount"],
                                "new biogenic CO2": exc["amount"],
                            }
                        )
                    else:
                        dataset["log parameters"].update(
                            {
                                "initial fossil CO2": exc["amount"],
                                "new fossil CO2": exc["amount"],
                            }
                        )

        return dataset

    def build_clinker_production_datasets(self) -> Dict[str, dict]:
        """
        Builds clinker production datasets for each IAM region.
        Adds CO2 capture and Storage, if needed.

        :return: a dictionary with IAM regions as keys and clinker production datasets as values.
        """

        # Fetch clinker production activities
        # and store them in a dictionary
        d_act_clinker = self.fetch_proxies(
            name="clinker production",
            ref_prod="clinker",
            production_variable="cement, dry feed rotary kiln",
        )

        for region, dataset in d_act_clinker.items():
            # calculate current thermal energy consumption per kg clinker
            energy_details = self.fetch_current_energy_details(dataset)
            current_energy_input_per_ton_clinker = sum(
                [d["energy"] for d in energy_details.values()]
            )

            # fetch the amount of biogenic CO2 emisisons
            bio_CO2 = sum(
                e["amount"]
                for e in dataset["exchanges"]
                if e["name"] == "Carbon dioxide, non-fossil"
            )

            # back-calculate the amount of waste fuel from
            # the biogenic CO2 emissions
            # biogenic CO2 / MJ for waste fuel
            waste_fuel_biogenic_co2_emission_factor = (
                self.fuels_specs["waste"]["co2"]
                * self.fuels_specs["waste"]["biogenic_share"]
            )

            waste_fuel_fossil_co2_emission_factor = self.fuels_specs["waste"]["co2"] * (
                1 - self.fuels_specs["waste"]["biogenic_share"]
            )
            # energy input of waste fuel in MJ
            energy_input_waste_fuel = bio_CO2 / waste_fuel_biogenic_co2_emission_factor
            # amount waste fuel, in kg
            amount_waste_fuel = (
                energy_input_waste_fuel / self.fuels_specs["waste"]["lhv"]
            )

            # add waste fuel to the energy details
            energy_details["market for waste plastic, mixture"] = {
                "amount": amount_waste_fuel,
                "energy": energy_input_waste_fuel * 1000,
                "fossil CO2": waste_fuel_fossil_co2_emission_factor
                * energy_input_waste_fuel,
                "biogenic CO2": bio_CO2,
            }

            # add the waste fuel energy input
            # to the total energy input
            current_energy_input_per_ton_clinker += energy_input_waste_fuel * 1000

            if "log parameters" not in dataset:
                dataset["log parameters"] = {}

            dataset["log parameters"].update(
                {
                    "initial energy input per ton clinker": current_energy_input_per_ton_clinker,
                }
            )

            # calculate the scaling factor
            # the correction factor applied to all fuel/electricity input is
            # equal to the ratio fuel/output in the year in question
            # divided by the ratio fuel/output in 2020

            scaling_factor = 1 / self.find_iam_efficiency_change(
                data=self.iam_data.cement_efficiencies,
                variable="cement, dry feed rotary kiln",
                location=dataset["location"],
            )

            # calculate new thermal energy
            # consumption per kg clinker
            new_energy_input_per_ton_clinker = (
                current_energy_input_per_ton_clinker * scaling_factor
            )

            # put a floor value of 3000 kj/kg clinker
            if new_energy_input_per_ton_clinker < 3000:
                new_energy_input_per_ton_clinker = 3000
            # and a ceiling value of 5000 kj/kg clinker
            elif new_energy_input_per_ton_clinker > 5000:
                new_energy_input_per_ton_clinker = 5000

            scaling_factor = (
                new_energy_input_per_ton_clinker / current_energy_input_per_ton_clinker
            )

            # rescale fuel consumption and emissions
            # rescale the fuel and electricity input
            dataset = self.rescale_fuel_inputs(dataset, scaling_factor, energy_details)

            # rescale combustion-related CO2 emissions
            dataset = self.rescale_emissions(dataset, energy_details, scaling_factor)

            # Carbon capture rate: share of capture of total CO2 emitted
            carbon_capture_rate = self.get_carbon_capture_rate(
                loc=dataset["location"], sector="cement, dry feed rotary kiln"
            )

            dataset["log parameters"].update(
                {
                    "carbon capture rate": carbon_capture_rate,
                }
            )

            # add CCS-related dataset
            if carbon_capture_rate > 0:
                # total CO2 emissions = bio CO2 emissions
                # + fossil CO2 emissions
                # + calcination emissions

                total_co2_emissions = dataset["log parameters"].get(
                    "new fossil CO2", 0
                ) + dataset["log parameters"].get("new biogenic CO2", 0)
                # share bio CO2 stored = sum of biogenic fuel emissions / total CO2 emissions
                bio_co2_stored = (
                    dataset["log parameters"].get("new biogenic CO2", 0)
                    / total_co2_emissions
                )

                # 0.11 kg CO2 leaks per kg captured
                # we need to align the CO2 composition with
                # the CO2 composition of the cement plant
                bio_co2_leaked = bio_co2_stored * 0.11

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
                    "name": "carbon dioxide, captured at cement production plant, with underground storage, post, 200 km",
                    "unit": "kilogram",
                    "location": dataset["location"],
                    "product": "carbon dioxide, captured and stored",
                }
                dataset["exchanges"].append(ccs_exc)

                # Update CO2 exchanges
                for exc in dataset["exchanges"]:
                    if (
                        exc["name"].lower().startswith("carbon dioxide")
                        and exc["type"] == "biosphere"
                    ):
                        exc["amount"] *= 1 - carbon_capture_rate

                        if "non-fossil" in exc["name"].lower():
                            dataset["log parameters"].update(
                                {
                                    "new biogenic CO2": exc["amount"],
                                }
                            )
                        else:
                            dataset["log parameters"].update(
                                {
                                    "new fossil CO2": exc["amount"],
                                }
                            )

            dataset["exchanges"] = [v for v in dataset["exchanges"] if v]

            # update comment
            dataset["comment"] = (
                "Dataset modified by `premise` based on IAM projections "
                + " for the cement industry.\n"
                + f"Calculated energy input per kg clinker: {np.round(new_energy_input_per_ton_clinker, 1) / 1000}"
                f" MJ/kg clinker.\n"
                + f"Rate of carbon capture: {int(carbon_capture_rate * 100)} pct.\n"
            ) + dataset["comment"]

        return d_act_clinker

    def add_datasets_to_database(self) -> None:
        """
        Runs a series of methods that create new clinker and cement production datasets
        and new cement market datasets.
        :return: Does not return anything. Modifies in place.
        """

        print("Start integration of cement data...")

        print("Create new clinker production datasets and delete old datasets")

        clinker_prod_datasets = list(self.build_clinker_production_datasets().values())
        self.database.extend(clinker_prod_datasets)

        # add to log
        for new_dataset in clinker_prod_datasets:
            self.write_log(new_dataset)
            # add it to list of created datasets
            self.modified_datasets[(self.model, self.scenario, self.year)][
                "created"
            ].append(
                (
                    new_dataset["name"],
                    new_dataset["reference product"],
                    new_dataset["location"],
                    new_dataset["unit"],
                )
            )

        print("Create new clinker market datasets and delete old datasets")
        clinker_market_datasets = list(
            self.fetch_proxies(
                name="market for clinker",
                ref_prod="clinker",
                production_variable="cement, dry feed rotary kiln",
            ).values()
        )

        self.database.extend(clinker_market_datasets)

        # add to log
        for new_dataset in clinker_market_datasets:
            self.write_log(new_dataset)
            # add it to list of created datasets
            self.modified_datasets[(self.model, self.scenario, self.year)][
                "created"
            ].append(
                (
                    new_dataset["name"],
                    new_dataset["reference product"],
                    new_dataset["location"],
                    new_dataset["unit"],
                )
            )

        print("Create new cement market datasets")

        # cement markets
        markets = ws.get_many(
            self.database,
            ws.contains("name", "market for cement"),
            ws.contains("reference product", "cement"),
            ws.doesnt_contain_any("name", ["factory", "tile", "sulphate", "plaster"]),
            ws.doesnt_contain_any("location", self.regions),
        )

        markets = {(m["name"], m["reference product"]) for m in markets}

        new_datasets = []

        for dataset in markets:
            new_cement_markets = self.fetch_proxies(
                name=dataset[0],
                ref_prod=dataset[1],
                production_variable="cement, dry feed rotary kiln",
            )

            # add to log
            for new_dataset in list(new_cement_markets.values()):
                self.write_log(new_dataset)
                # add it to list of created datasets
                self.modified_datasets[(self.model, self.scenario, self.year)][
                    "created"
                ].append(
                    (
                        new_dataset["name"],
                        new_dataset["reference product"],
                        new_dataset["location"],
                        new_dataset["unit"],
                    )
                )

            new_datasets.extend(list(new_cement_markets.values()))

        self.database.extend(new_datasets)

        print(
            "Create new cement production datasets and "
            "adjust electricity consumption"
        )
        # cement production
        production = ws.get_many(
            self.database,
            ws.contains("name", "cement production"),
            ws.contains("reference product", "cement"),
            ws.doesnt_contain_any("name", ["factory", "tile", "sulphate", "plaster"]),
        )

        production = {(p["name"], p["reference product"]) for p in production}

        new_datasets = []

        for dataset in production:
            # Fetch proxy datasets (one per IAM region)
            # Delete old datasets
            new_cement_production = self.fetch_proxies(
                name=dataset[0],
                ref_prod=dataset[1],
                production_variable="cement, dry feed rotary kiln",
            )

            # add to log
            for new_dataset in list(new_cement_production.values()):
                self.write_log(dataset=new_dataset, status="updated")
                # add it to list of created datasets
                self.modified_datasets[(self.model, self.scenario, self.year)][
                    "created"
                ].append(
                    (
                        new_dataset["name"],
                        new_dataset["reference product"],
                        new_dataset["location"],
                        new_dataset["unit"],
                    )
                )

            new_datasets.extend(list(new_cement_production.values()))

        self.database.extend(new_datasets)

        print("Done!")

    def write_log(self, dataset, status="created"):
        """
        Write log file.
        """

        logger.info(
            f"{status}|{self.model}|{self.scenario}|{self.year}|"
            f"{dataset['name']}|{dataset['location']}|"
            f"{dataset.get('log parameters', {}).get('initial energy input per ton clinker', '')}|"
            f"{dataset.get('log parameters', {}).get('new energy input per ton clinker', '')}|"
            f"{dataset.get('log parameters', {}).get('carbon capture rate', '')}|"
            f"{dataset.get('log parameters', {}).get('initial fossil CO2', '')}|"
            f"{dataset.get('log parameters', {}).get('initial biogenic CO2', '')}|"
            f"{dataset.get('log parameters', {}).get('new fossil CO2', '')}|"
            f"{dataset.get('log parameters', {}).get('new biogenic CO2', '')}|"
            f"{dataset.get('log parameters', {}).get('electricity generated', '')}|"
            f"{dataset.get('log parameters', {}).get('electricity consumed', '')}|"
            f"{dataset.get('log parameters', {}).get('old clinker-to-cement ratio', '')}|"
            f"{dataset.get('log parameters', {}).get('new clinker-to-cement ratio', '')}"
        )
