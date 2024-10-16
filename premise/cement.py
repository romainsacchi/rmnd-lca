"""
cement.py contains the class `Cement`, which inherits from `BaseTransformation`.
This class transforms the cement markets and clinker and cement production activities
of the wurst database, based on projections from the IAM scenario.
It eventually re-links all the cement-consuming activities (e.g., concrete production)
of the wurst database to the newly created cement markets.

"""

import copy
import uuid

from .export import biosphere_flows_dictionary
from .logger import create_logger
from .transformation import (
    BaseTransformation,
    IAMDataCollection,
    InventorySet,
    List,
    np,
    ws,
)
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
        cement.replace_clinker_production_with_markets()
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

        mapping = InventorySet(self.database)
        self.cement_fuels_map: dict = mapping.generate_cement_fuels_map()

        # reverse the fuel map to get a mapping from ecoinvent to premise
        self.fuel_map_reverse: dict = {}

        self.fuel_map: dict = mapping.generate_fuel_map()

        for key, value in self.fuel_map.items():
            for v in list(value):
                self.fuel_map_reverse[v] = key

        self.biosphere_dict = biosphere_flows_dictionary(self.version)

    def build_CCS_datasets(self):
        ccs_datasets = {
            "cement, dry feed rotary kiln, efficient, with on-site CCS": {
                "name": "carbon dioxide, captured at cement production plant, using direct separation",
                "reference product": "carbon dioxide, captured at cement plant",
            },
            "cement, dry feed rotary kiln, efficient, with oxyfuel CCS": {
                "name": "carbon dioxide, captured at cement production plant, using oxyfuel",
                "reference product": "carbon dioxide, captured at cement plant",
            },
            "cement, dry feed rotary kiln, efficient, with MEA CCS": {
                "name": "carbon dioxide, captured at cement production plant, using monoethanolamine",
                "reference product": "carbon dioxide, captured at cement plant",
            },
        }

        for variable in ccs_datasets:
            datasets = self.fetch_proxies(
                name=ccs_datasets[variable]["name"],
                ref_prod=ccs_datasets[variable]["reference product"],
            )

            if variable == "cement, dry feed rotary kiln, efficient, with MEA CCS":
                # we adjust the heat needs by subtraction 3.66 MJ with what
                # the plant is expected to produce as excess heat

                # Heat, as steam: 3.66 MJ/kg CO2 captured in 2020,
                # decreasing to 2.6 GJ/t by 2050, by looking at
                # the best-performing state-of-the-art technologies today
                # https://www.globalccsinstitute.com/wp-content/uploads/2022/05/State-of-the-Art-CCS-Technologies-2022.pdf
                # minus excess heat generated on site
                # the contribution of excess heat is assumed to be
                # 30% of heat requirement.

                heat_input = np.clip(
                    np.interp(self.year, [2020, 2050], [3.66, 2.6]), 2.6, 3.66
                )
                excess_heat_generation = 0.3  # 30%
                fossil_heat_input = heat_input - (excess_heat_generation * heat_input)

                for region, dataset in datasets.items():
                    for exc in ws.technosphere(
                        dataset, ws.contains("unit", "megajoule")
                    ):
                        exc["amount"] = fossil_heat_input

            for dataset in datasets.values():
                self.add_to_index(dataset)
                self.write_log(dataset)
                self.database.append(dataset)

        # also create region-specific air separation datasets
        datasets_to_regionalize = [
            (
                "industrial gases production, cryogenic air separation"
                if self.version == "3.10"
                else "air separation, cryogenic"
            ),
            "market for oxygen, liquid",
        ]

        for ds_to_regionlaize in datasets_to_regionalize:

            air_separation = self.fetch_proxies(
                name=ds_to_regionlaize,
                ref_prod="oxygen, liquid",
            )

            for dataset in air_separation.values():
                self.add_to_index(dataset)
                self.write_log(dataset)
                self.database.append(dataset)

    def build_clinker_production_datasets(self) -> list:
        """
        Builds clinker production datasets for each IAM region.
        Adds CO2 capture and Storage, if needed.

        :return: a dictionary with IAM regions as keys
        and clinker production datasets as values.
        """

        variables = [
            "cement, dry feed rotary kiln",
            "cement, dry feed rotary kiln, efficient",
            "cement, dry feed rotary kiln, efficient, with on-site CCS",
            "cement, dry feed rotary kiln, efficient, with oxyfuel CCS",
            "cement, dry feed rotary kiln, efficient, with MEA CCS",
        ]

        datasets = []

        # Fetch clinker production activities
        # and store them in a dictionary
        clinker = self.fetch_proxies(
            name="clinker production",
            ref_prod="clinker",
            production_variable="cement, dry feed rotary kiln",
            geo_mapping={r: "Europe without Switzerland" for r in self.regions},
        )

        for variable in variables:
            if variable in self.iam_data.cement_markets.coords["variables"].values:

                d_act_clinker = copy.deepcopy(clinker)
                # remove `code` field
                for region, dataset in d_act_clinker.items():
                    dataset["code"] = uuid.uuid4().hex
                    for exc in ws.production(dataset):
                        if "input" in exc:
                            del exc["input"]

                if variable != "cement, dry feed rotary kiln":
                    # rename datasets
                    for region, dataset in d_act_clinker.items():
                        dataset["name"] = (
                            f"{dataset['name']}, {variable.replace('cement, dry feed rotary kiln, ', '')}"
                        )
                        for e in dataset["exchanges"]:
                            if e["type"] == "production":
                                e["name"] = (
                                    f"{e['name']}, {variable.replace('cement, dry feed rotary kiln, ', '')}"
                                )

                for region, dataset in d_act_clinker.items():

                    # from Kellenberger at al. 2007, the total energy
                    # input per ton of clinker is 3.4 GJ/ton clinker
                    current_energy_input_per_ton_clinker = 3400

                    # calculate the scaling factor
                    # the correction factor applied to hard coal input
                    # we assume that any fuel use reduction would in priority
                    # affect hard coal use

                    scaling_factor = 1 / self.find_iam_efficiency_change(
                        data=self.iam_data.cement_efficiencies,
                        variable=variable,
                        location=dataset["location"],
                    )

                    new_energy_input_per_ton_clinker = 3400

                    if "log parameters" not in dataset:
                        dataset["log parameters"] = {}

                    dataset["log parameters"][
                        "initial energy input per ton clinker"
                    ] = current_energy_input_per_ton_clinker

                    if not np.isnan(scaling_factor):
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

                        # but if efficient kiln,
                        # set the energy input to 3000 kJ/kg clinker
                        if variable.startswith(
                            "cement, dry feed rotary kiln, efficient"
                        ):
                            new_energy_input_per_ton_clinker = 3000

                        dataset["log parameters"][
                            "new energy input per ton clinker"
                        ] = int(new_energy_input_per_ton_clinker)

                        scaling_factor = (
                            new_energy_input_per_ton_clinker
                            / current_energy_input_per_ton_clinker
                        )

                        dataset["log parameters"][
                            "energy scaling factor"
                        ] = scaling_factor

                        # rescale hard coal consumption and related emissions
                        coal_specs = self.fuels_specs["hard coal"]
                        old_coal_input, new_coal_input = 0, 0
                        for exc in ws.technosphere(
                            dataset,
                            ws.contains("name", "hard coal"),
                        ):
                            # in kJ
                            old_coal_input = float(exc["amount"] * coal_specs["lhv"])
                            # in MJ
                            new_coal_input = old_coal_input - (
                                (
                                    current_energy_input_per_ton_clinker
                                    - new_energy_input_per_ton_clinker
                                )
                                / 1000
                            )
                            exc["amount"] = np.clip(
                                new_coal_input / coal_specs["lhv"], 0, None
                            )

                        # rescale combustion-related fossil CO2 emissions
                        for exc in ws.biosphere(
                            dataset,
                            ws.contains("name", "Carbon dioxide"),
                        ):
                            if exc["name"] == "Carbon dioxide, fossil":
                                dataset["log parameters"]["initial fossil CO2"] = exc[
                                    "amount"
                                ]
                                co2_reduction = (
                                    old_coal_input - new_coal_input
                                ) * coal_specs["co2"]
                                exc["amount"] -= co2_reduction
                                dataset["log parameters"]["new fossil CO2"] = exc[
                                    "amount"
                                ]

                            if exc["name"] == "Carbon dioxide, non-fossil":
                                dataset["log parameters"]["initial biogenic CO2"] = exc[
                                    "amount"
                                ]

                    # add 0.005 kg/kg clinker of ammonia use for NOx removal
                    # according to Muller et al., 2024
                    for exc in ws.technosphere(
                        dataset,
                        ws.contains("name", "market for ammonia"),
                    ):
                        if (
                            variable
                            == "cement, dry feed rotary kiln, efficient, with MEA CCS"
                        ):
                            exc["amount"] = 0.00662
                        else:
                            exc["amount"] = 0.005

                    # reduce NOx emissions
                    # according to Muller et al., 2024
                    for exc in ws.biosphere(
                        dataset,
                        ws.contains("name", "Nitrogen oxides"),
                    ):
                        if variable in [
                            "cement, dry feed rotary kiln, efficient, with on-site CCS",
                            "cement, dry feed rotary kiln, efficient, with oxyfuel CCS",
                        ]:
                            exc["amount"] = 1.22e-5
                        elif (
                            variable
                            == "cement, dry feed rotary kiln, efficient, with MEA CCS"
                        ):
                            exc["amount"] = 3.8e-4
                        else:
                            exc["amount"] = 7.6e-4

                    # reduce Mercury and SOx emissions
                    # according to Muller et al., 2024
                    if variable in [
                        "cement, dry feed rotary kiln, efficient, with on-site CCS",
                        "cement, dry feed rotary kiln, efficient, with oxyfuel CCS",
                        "cement, dry feed rotary kiln, efficient, with MEA CCS",
                    ]:
                        for exc in ws.biosphere(
                            dataset,
                            ws.either(
                                *[
                                    ws.contains("name", name)
                                    for name in [
                                        "Mercury",
                                        "Sulfur dioxide",
                                    ]
                                ]
                            ),
                        ):
                            exc["amount"] *= 1 - 0.999

                    if self.model == "image":
                        # add CCS datasets
                        ccs_datasets = {
                            "cement, dry feed rotary kiln, efficient, with on-site CCS": {
                                "name": "carbon dioxide, captured at cement production plant, using direct separation",
                                "reference product": "carbon dioxide, captured at cement plant",
                                "capture share": 0.95,  # 95% of process emissions (calcination) are captured
                            },
                            "cement, dry feed rotary kiln, efficient, with oxyfuel CCS": {
                                "name": "carbon dioxide, captured at cement production plant, using oxyfuel",
                                "reference product": "carbon dioxide, captured at cement plant",
                                "capture share": 0.9,
                            },
                            "cement, dry feed rotary kiln, efficient, with MEA CCS": {
                                "name": "carbon dioxide, captured at cement production plant, using monoethanolamine",
                                "reference product": "carbon dioxide, captured at cement plant",
                                "capture share": 0.9,
                            },
                        }

                        if variable in ccs_datasets:
                            CO2_amount = sum(
                                e["amount"]
                                for e in ws.biosphere(
                                    dataset,
                                    ws.contains("name", "Carbon dioxide"),
                                )
                            )
                            if (
                                variable
                                == "cement, dry feed rotary kiln, efficient, with on-site CCS"
                            ):
                                # only 95% of process emissions (calcination) are captured
                                CCS_amount = (
                                    0.543 * ccs_datasets[variable]["capture share"]
                                )
                            else:
                                CCS_amount = (
                                    CO2_amount * ccs_datasets[variable]["capture share"]
                                )

                            dataset["log parameters"]["carbon capture rate"] = (
                                CCS_amount / CO2_amount
                            )

                            ccs_exc = {
                                "uncertainty type": 0,
                                "loc": CCS_amount,
                                "amount": CCS_amount,
                                "type": "technosphere",
                                "production volume": 0,
                                "name": ccs_datasets[variable]["name"],
                                "unit": "kilogram",
                                "location": dataset["location"],
                                "product": ccs_datasets[variable]["reference product"],
                            }
                            dataset["exchanges"].append(ccs_exc)

                            # Update CO2 exchanges
                            for exc in ws.biosphere(
                                dataset,
                                ws.contains("name", "Carbon dioxide, fossil"),
                            ):
                                if (
                                    variable
                                    != "cement, dry feed rotary kiln, efficient, with on-site CCS"
                                ):
                                    exc["amount"] *= (
                                        CO2_amount - CCS_amount
                                    ) / CO2_amount
                                else:
                                    exc["amount"] -= CCS_amount

                                # make sure it's not negative
                                if exc["amount"] < 0:
                                    exc["amount"] = 0

                                dataset["log parameters"]["new fossil CO2"] = exc[
                                    "amount"
                                ]

                            # Update biogenic CO2 exchanges
                            if (
                                variable
                                != "cement, dry feed rotary kiln, efficient, with on-site CCS"
                            ):
                                for exc in ws.biosphere(
                                    dataset,
                                    ws.contains("name", "Carbon dioxide, non-fossil"),
                                ):
                                    dataset["log parameters"][
                                        "initial biogenic CO2"
                                    ] = exc["amount"]
                                    exc["amount"] *= (
                                        CO2_amount - CCS_amount
                                    ) / CO2_amount

                                    # make sure it's not negative
                                    if exc["amount"] < 0:
                                        exc["amount"] = 0

                                    dataset["log parameters"]["new biogenic CO2"] = exc[
                                        "amount"
                                    ]

                                    biogenic_CO2_reduction = (
                                        dataset["log parameters"][
                                            "initial biogenic CO2"
                                        ]
                                        - dataset["log parameters"]["new biogenic CO2"]
                                    )
                                    # add a flow of "Carbon dioxide, in air" to reflect
                                    # the permanent storage of biogenic CO2
                                    dataset["exchanges"].append(
                                        {
                                            "uncertainty type": 0,
                                            "loc": biogenic_CO2_reduction,
                                            "amount": biogenic_CO2_reduction,
                                            "type": "biosphere",
                                            "name": "Carbon dioxide, in air",
                                            "unit": "kilogram",
                                            "categories": (
                                                "natural resource",
                                                "in air",
                                            ),
                                            "comment": "Permanent storage of biogenic CO2",
                                            "input": (
                                                "biosphere3",
                                                self.biosphere_dict[
                                                    (
                                                        "Carbon dioxide, in air",
                                                        "natural resource",
                                                        "in air",
                                                        "kilogram",
                                                    )
                                                ],
                                            ),
                                        }
                                    )

                    else:
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
                        if (
                            not np.isnan(carbon_capture_rate)
                            and carbon_capture_rate > 0
                        ):
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

                            # add an input from this CCS dataset in the clinker dataset
                            ccs_exc = {
                                "uncertainty type": 0,
                                "loc": 0,
                                "amount": float(
                                    total_co2_emissions * carbon_capture_rate
                                ),
                                "type": "technosphere",
                                "production volume": 0,
                                "name": "carbon dioxide, captured at cement production plant, using monoethanolamine",
                                "unit": "kilogram",
                                "location": dataset["location"],
                                "product": "carbon dioxide, captured at cement plant",
                            }

                            # add an input from the CCS dataset in the clinker dataset
                            # and add it to the database
                            dataset["exchanges"].append(ccs_exc)

                            # Update CO2 exchanges
                            for exc in ws.biosphere(
                                dataset,
                                ws.contains("name", "Carbon dioxide"),
                            ):
                                if exc["name"] == "Carbon dioxide, fossil":
                                    exc["amount"] *= 1 - carbon_capture_rate
                                    dataset["log parameters"]["new fossil CO2"] = exc[
                                        "amount"
                                    ]

                                if exc["name"] == "Carbon dioxide, non-fossil":
                                    exc["amount"] *= 1 - carbon_capture_rate
                                    dataset["log parameters"]["new biogenic CO2"] = exc[
                                        "amount"
                                    ]

                            # add a flow of "Carbon dioxide, in air" to reflect
                            # the permanent storage of biogenic CO2
                            biogenic_CO2_reduction = dataset["log parameters"].get(
                                "initial biogenic CO2", 0.0
                            ) - dataset["log parameters"].get("new biogenic CO2", 0.0)
                            dataset["exchanges"].append(
                                {
                                    "uncertainty type": 0,
                                    "loc": biogenic_CO2_reduction,
                                    "amount": biogenic_CO2_reduction,
                                    "type": "biosphere",
                                    "name": "Carbon dioxide, in air",
                                    "unit": "kilogram",
                                    "categories": (
                                        "natural resource",
                                        "in air",
                                    ),
                                    "comment": "Permanent storage of biogenic CO2",
                                    "input": (
                                        "biosphere3",
                                        self.biosphere_dict[
                                            (
                                                "Carbon dioxide, in air",
                                                "natural resource",
                                                "in air",
                                                "kilogram",
                                            )
                                        ],
                                    ),
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

                    datasets.append(dataset)

        return datasets

    def replace_clinker_production_with_markets(self):
        """
        Some cement production datasets in ecoinvent receive an input from clinker production datasets.
        This is problematic because it will not benefit from the new cement markets, containing alternative clinker production pathways.
        So we replace the clinker production datasets with the clinker markets.
        """

        for ds in ws.get_many(
            self.database,
            ws.contains("name", "cement production"),
            ws.contains("reference product", "cement"),
            ws.equals("unit", "kilogram"),
        ):
            for exc in ws.technosphere(ds):
                if exc["name"] == "clinker production" and exc["product"] == "clinker":
                    exc["name"] = "market for clinker"

    def add_datasets_to_database(self) -> None:
        """
        Runs a series of methods that create new clinker and cement production datasets
        and new cement market datasets.
        :return: Does not return anything. Modifies in place.
        """

        # create CCS datasets
        self.build_CCS_datasets()
        clinker_prod_datasets = self.build_clinker_production_datasets()
        self.database.extend(clinker_prod_datasets)

        # add to log
        for new_dataset in clinker_prod_datasets:
            self.write_log(new_dataset)
            # add it to list of created datasets
            self.add_to_index(new_dataset)

        variables = [
            "cement, dry feed rotary kiln",
            "cement, dry feed rotary kiln, efficient",
            "cement, dry feed rotary kiln, efficient, with on-site CCS",
            "cement, dry feed rotary kiln, efficient, with oxyfuel CCS",
            "cement, dry feed rotary kiln, efficient, with MEA CCS",
        ]

        clinker_market_datasets = self.fetch_proxies(
            name="market for clinker",
            ref_prod="clinker",
            production_variable=[
                v
                for v in variables
                if v in self.iam_data.cement_markets.coords["variables"].values
            ],
        )

        clinker_market_datasets = {
            k: v
            for k, v in clinker_market_datasets.items()
            if self.iam_data.cement_markets.sel(region=k)
            .sum(dim="variables")
            .interp(year=self.year)
            > 0
        }

        for region, ds in clinker_market_datasets.items():
            ds["exchanges"] = [
                v
                for v in ds["exchanges"]
                if v["type"] == "production" or v["unit"] == "ton kilometer"
            ]
            for variable in variables:
                if variable in self.iam_data.cement_markets.coords["variables"].values:
                    if self.year in self.iam_data.cement_markets.coords["year"].values:
                        share = self.iam_data.cement_markets.sel(
                            variables=variable, region=region, year=self.year
                        ).values
                    else:
                        share = (
                            self.iam_data.cement_markets.sel(
                                variables=variable, region=region
                            )
                            .interp(year=self.year)
                            .values
                        )

                    if share > 0:
                        if variable == "cement, dry feed rotary kiln":
                            name = "clinker production"
                        else:
                            name = f"clinker production, {variable.replace('cement, dry feed rotary kiln, ', '')}"
                        ds["exchanges"].append(
                            {
                                "uncertainty type": 0,
                                "loc": float(share),
                                "amount": float(share),
                                "type": "technosphere",
                                "production volume": 0,
                                "name": name,
                                "unit": "kilogram",
                                "location": region,
                                "product": "clinker",
                            }
                        )

        self.database.extend(clinker_market_datasets.values())

        # add to log
        for new_dataset in clinker_market_datasets.values():
            self.write_log(new_dataset)
            # add it to list of created datasets
            self.add_to_index(new_dataset)

        # exclude the regionalization of these datasets
        # because they are very rarely used in the database
        excluded = [
            "factory",
            "tile",
            "sulphate",
            "plaster",
            "Portland Slag",
            "CP II-Z",
            "CP IV",
            "CP V RS",
            "Portland SR3",
            "CEM II/A-S",
            "CEM II/A-V",
            "CEM II/B-L",
            "CEM II/B-S",
            "type I (SM)",
            "type I-PM",
            "type IP/P",
            "type IS",
            "type S",
            "CEM III/C",
            "CEM V/A",
            "CEM V/B",
            "CEM II/A-L",
            "CEM III/B",
            "Pozzolana Portland",
            "ART",
            "type IP",
            "CEM IV/A",
            "CEM IV/B",
            "type ICo",
            "carbon",
        ]

        # cement markets
        markets = list(
            ws.get_many(
                self.database,
                ws.contains("name", "market for cement"),
                ws.contains("reference product", "cement"),
                ws.doesnt_contain_any(
                    "name",
                    excluded,
                ),
                ws.doesnt_contain_any("location", self.regions),
            )
        )

        unique_markets = list(
            set([(m["name"], m["reference product"]) for m in markets])
        )

        new_datasets = []

        for dataset in unique_markets:
            new_cement_markets = self.fetch_proxies(
                name=dataset[0],
                ref_prod=dataset[1],
                production_variable="cement, dry feed rotary kiln",
                subset=markets,
            )

            # add to log
            for new_dataset in list(new_cement_markets.values()):
                self.write_log(new_dataset)
                # add it to list of created datasets
                self.add_to_index(new_dataset)

            new_datasets.extend(list(new_cement_markets.values()))

        self.database.extend(new_datasets)

        # cement production
        production = list(
            ws.get_many(
                self.database,
                ws.contains("name", "cement production"),
                ws.contains("reference product", "cement"),
                ws.doesnt_contain_any("name", excluded),
            )
        )

        reduced_production = list(
            set([(p["name"], p["reference product"]) for p in production])
        )

        new_datasets = []

        for dataset in reduced_production:
            # Fetch proxy datasets (one per IAM region)
            # Delete old datasets
            new_cement_production = self.fetch_proxies(
                name=dataset[0],
                ref_prod=dataset[1],
                production_variable="cement, dry feed rotary kiln",
                subset=production,
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
            f"{dataset.get('log parameters', {}).get('energy scaling factor', '')}|"
            f"{dataset.get('log parameters', {}).get('new energy input per ton clinker', '')}|"
            f"{dataset.get('log parameters', {}).get('carbon capture rate', '')}|"
            f"{dataset.get('log parameters', {}).get('initial fossil CO2', '')}|"
            f"{dataset.get('log parameters', {}).get('initial biogenic CO2', '')}|"
            f"{dataset.get('log parameters', {}).get('new fossil CO2', '')}|"
            f"{dataset.get('log parameters', {}).get('new biogenic CO2', '')}|"
            f"{dataset.get('log parameters', {}).get('electricity generated', '')}|"
            f"{dataset.get('log parameters', {}).get('electricity consumed', '')}"
        )
