"""
cement.py contains the class `Cement`, which inherits from `BaseTransformation`.
This class transforms the cement markets and clinker and cement production activities
of the wurst database, based on projections from the IAM scenario.
It eventually re-links all the cement-consuming activities (e.g., concrete production)
of the wurst database to the newly created cement markets.

"""

from collections import defaultdict

from .logger import create_logger
from .transformation import BaseTransformation, Dict, IAMDataCollection, List, np, ws
from .validation import CementValidation

logger = create_logger("cement")


def _update_cement(scenario, version, system_model):
    cement = Cement(
        database=scenario["database"],
        model=scenario["model"],
        pathway=scenario["pathway"],
        iam_data=scenario["iam data"],
        year=scenario["year"],
        version=version,
        system_model=system_model,
        cache=scenario.get("cache"),
        index=scenario.get("index"),
    )

    if scenario["iam data"].cement_markets is not None:
        cement.add_datasets_to_database()
        cement.relink_datasets()
        scenario["database"] = cement.database
        scenario["index"] = cement.index
        scenario["cache"] = cement.cache

        validate = CementValidation(
            model=scenario["model"],
            scenario=scenario["pathway"],
            year=scenario["year"],
            regions=scenario["iam data"].regions,
            database=cement.database,
            iam_data=scenario["iam data"],
        )

        validate.run_cement_checks()
    else:
        print("No cement markets found in IAM data. Skipping.")

    return scenario


class Cement(BaseTransformation):
    """
    Class that modifies clinker and cement production datasets in ecoinvent.
    It creates region-specific new clinker production datasets (and deletes the original ones).
    It adjusts the kiln efficiency based on the improvement indicated in the IAM file, relative to 2020.
    It adds CCS, if indicated in the IAM file.
    It creates regions-specific cement production datasets (and deletes the original ones).
    It adjusts electricity consumption in cement production datasets.
    It creates regions-specific cement market datasets (and deletes the original ones).


    :ivar database: wurst database, which is a list of dictionaries
    :ivar iam_data: IAM data
    :ivar model: name of the IAM model (e.g., "remind", "image")
    :ivar pathway: name of the IAM scenario (e.g., "SSP2-19")
    :ivar year: year of the pathway (e.g., 2030)
    :ivar version: version of ecoinvent database (e.g., "3.7")
    :ivar system_model: name of the system model (e.g., "attributional", "consequential")

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
        self.version = version

    def fetch_current_energy_details(self, dataset):
        # Fetches the current energy consumption per ton of clinker

        d_fuels = defaultdict(
            lambda: {"amount": 0, "energy": 0, "fossil CO2": 0, "biogenic CO2": 0}
        )

        for exc in ws.technosphere(dataset):
            fuel_name = exc["name"]

            if fuel_name in self.cement_fuels_map["cement, dry feed rotary kiln"]:
                fuel_data = self.fuels_specs[self.fuel_map_reverse[fuel_name]]
                co2_emission = fuel_data["co2"]
                biogenic_share = fuel_data["biogenic_share"]

                # Calculate the energy once for the given exc
                input_energy = (
                    self.calculate_input_energy(exc["name"], exc["amount"], exc["unit"])
                    * 1000
                )

                # Update the dictionary in one go
                d_fuels[fuel_name]["amount"] += exc["amount"]
                d_fuels[fuel_name]["energy"] += input_energy
                d_fuels[fuel_name]["fossil CO2"] += (
                    co2_emission * input_energy * (1 - biogenic_share)
                )
                d_fuels[fuel_name]["biogenic CO2"] += (
                    co2_emission * input_energy * biogenic_share
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
            e["amount"]
            for e in dataset["exchanges"]
            if e["name"] == "Carbon dioxide, non-fossil"
        )
        non_bio_co2 = sum(
            e["amount"]
            for e in dataset["exchanges"]
            if e["name"] == "Carbon dioxide, fossil"
        )

        return bio_co2 / (bio_co2 + non_bio_co2)

    def rescale_fuel_inputs(self, dataset, scaling_factor, energy_details):
        if scaling_factor != 1 and scaling_factor > 0:
            for exc in dataset["exchanges"]:
                if exc["name"] in self.cement_fuels_map["cement, dry feed rotary kiln"]:
                    exc["amount"] *= scaling_factor

                    # update energy_details
                    energy_details[exc["name"]]["amount"] *= scaling_factor
                    energy_details[exc["name"]]["energy"] *= scaling_factor
                    energy_details[exc["name"]]["fossil CO2"] *= scaling_factor
                    energy_details[exc["name"]]["biogenic CO2"] *= scaling_factor

        new_energy_input = sum(d["energy"] for d in energy_details.values())

        dataset["log parameters"].update(
            {
                "new energy input per ton clinker": new_energy_input,
            }
        )

        return dataset

    def rescale_emissions(self, dataset: dict, energy_details: dict) -> dict:
        for exc in ws.biosphere(dataset, ws.contains("name", "Carbon dioxide")):
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

                exc["amount"] = total_biogenic_CO2 / 1000
                dataset["log parameters"].update(
                    {
                        "new biogenic CO2": total_biogenic_CO2 / 1000,
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
                exc["amount"] = 0.525 + (total_fossil_CO2 / 1000)
                dataset["log parameters"].update(
                    {
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
                d["energy"] for d in energy_details.values()
            )

            # fetch the amount of biogenic CO2 emissions
            bio_CO2 = sum(
                e["amount"]
                for e in ws.biosphere(
                    dataset, ws.contains("name", "Carbon dioxide, non-fossil")
                )
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

            # add the waste fuel input to the dataset
            if amount_waste_fuel != 0:
                dataset["exchanges"].append(
                    {
                        "uncertainty type": 0,
                        "loc": 0,
                        "amount": amount_waste_fuel * -1,
                        "type": "technosphere",
                        "production volume": 0,
                        "name": "clinker production",
                        "unit": "kilogram",
                        "location": "RoW",
                        "product": "waste plastic, mixture",
                    }
                )

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

            new_energy_input_per_ton_clinker = 0

            if not np.isnan(scaling_factor) and scaling_factor > 0.0:
                # calculate new thermal energy
                # consumption per kg clinker
                new_energy_input_per_ton_clinker = (
                    current_energy_input_per_ton_clinker * scaling_factor
                )

                # put a floor value of 3100 kj/kg clinker
                if new_energy_input_per_ton_clinker < 3100:
                    new_energy_input_per_ton_clinker = 3100
                # and a ceiling value of 5000 kj/kg clinker
                elif new_energy_input_per_ton_clinker > 5000:
                    new_energy_input_per_ton_clinker = 5000

                scaling_factor = (
                    new_energy_input_per_ton_clinker
                    / current_energy_input_per_ton_clinker
                )

                # rescale fuel consumption and emissions
                # rescale the fuel and electricity input
                dataset = self.rescale_fuel_inputs(
                    dataset, scaling_factor, energy_details
                )

                # rescale combustion-related CO2 emissions
                dataset = self.rescale_emissions(dataset, energy_details)

            # Carbon capture rate: share of capture of total CO2 emitted
            carbon_capture_rate = self.get_carbon_capture_rate(
                loc=dataset["location"], sector="cement"
            )

            # add 10% loss
            carbon_capture_rate *= 0.9

            dataset["log parameters"].update(
                {
                    "carbon capture rate": float(carbon_capture_rate),
                }
            )

            # add CCS-related dataset
            if not np.isnan(carbon_capture_rate) and carbon_capture_rate > 0:
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
                    "amount": float(total_co2_emissions * carbon_capture_rate),
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

        clinker_prod_datasets = list(self.build_clinker_production_datasets().values())
        self.database.extend(clinker_prod_datasets)

        # add to log
        for new_dataset in clinker_prod_datasets:
            self.write_log(new_dataset)
            # add it to list of created datasets
            self.add_to_index(new_dataset)

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
            self.add_to_index(new_dataset)

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
                self.add_to_index(new_dataset)

            new_datasets.extend(list(new_cement_markets.values()))

        self.database.extend(new_datasets)

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
                self.add_to_index(new_dataset)

            new_datasets.extend(list(new_cement_production.values()))

        self.database.extend(new_datasets)

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
            f"{dataset.get('log parameters', {}).get('electricity consumed', '')}"
        )
