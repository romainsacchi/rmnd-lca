"""
electricity.py contains the class `Electricity`, which inherits from `BaseTransformation`.
This class transforms the electricity markets and power plants of the wurst database,
based on projections from the IAM scenario.
It also creates electricity markets which mix is a weighted-average
over a certain period (e.g., 10, 20 years).
It eventually re-links all the electricity-consuming activities of the wurst database to
the newly created electricity markets.

"""

import logging.config
import re
from collections import defaultdict
from functools import lru_cache
from pathlib import Path

import wurst
import yaml

from . import VARIABLES_DIR
from .data_collection import get_delimiter
from .export import biosphere_flows_dictionary
from .transformation import (
    BaseTransformation,
    Dict,
    IAMDataCollection,
    InventorySet,
    List,
    Tuple,
    get_shares_from_production_volume,
    get_suppliers_of_a_region,
    np,
    uuid,
    ws,
)
from .utils import DATA_DIR, eidb_label, get_efficiency_solar_photovoltaics

LOSS_PER_COUNTRY = DATA_DIR / "electricity" / "losses_per_country.csv"
IAM_BIOMASS_VARS = VARIABLES_DIR / "biomass_variables.yaml"

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

logger = logging.getLogger("electricity")


def get_losses_per_country_dict() -> Dict[str, Dict[str, float]]:
    """
    Create a dictionary with ISO country codes as keys and loss ratios as values.
    :return: ISO country code (keys) to loss ratio (values) dictionary
    :rtype: dict
    """

    if not LOSS_PER_COUNTRY.is_file():
        raise FileNotFoundError(
            "The production per country dictionary " "file could not be found."
        )

    sep = get_delimiter(filepath=LOSS_PER_COUNTRY)

    with open(LOSS_PER_COUNTRY, encoding="utf-8") as file:
        csv_list = [[val.strip() for val in r.split(sep)] for r in file.readlines()]

    (_, *header), *data = csv_list
    csv_dict = {}
    for row in data:
        key, *values = row
        csv_dict[key] = {key: float(value) for key, value in zip(header, values)}

    return csv_dict


def get_production_weighted_losses(
    losses: Dict[str, Dict[str, float]], locs: List[str]
) -> Dict[str, Dict[str, float]]:
    """
    Return the transformation, transmission and distribution losses at a given voltage level for a given location.
    A weighted average is made of the locations contained in the IAM region.

    """

    # Fetch locations contained in IAM region

    cumul_prod, transf_loss = 0.0, 0.0
    for loc in locs:
        dict_loss = losses.get(
            loc,
            {"Transformation loss high voltage": 0.0, "Production volume": 0.0},
        )

        transf_loss += (
            dict_loss["Transformation loss high voltage"]
            * dict_loss["Production volume"]
        )
        cumul_prod += dict_loss["Production volume"]
    transf_loss /= cumul_prod

    high = {"transf_loss": transf_loss, "distr_loss": 0.0}

    cumul_prod, transf_loss, distr_loss = 0.0, 0.0, 0.0

    for loc in locs:
        dict_loss = losses.get(
            loc,
            {
                "Transformation loss medium voltage": 0,
                "Transmission loss to medium voltage": 0,
                "Production volume": 0,
            },
        )
        transf_loss += (
            dict_loss["Transformation loss medium voltage"]
            * dict_loss["Production volume"]
        )
        distr_loss += (
            dict_loss["Transmission loss to medium voltage"]
            * dict_loss["Production volume"]
        )
        cumul_prod += dict_loss["Production volume"]
    transf_loss /= cumul_prod
    distr_loss /= cumul_prod

    medium = {"transf_loss": transf_loss, "distr_loss": distr_loss}

    cumul_prod, transf_loss, distr_loss = 0.0, 0.0, 0.0

    for loc in locs:
        dict_loss = losses.get(
            loc,
            {
                "Transformation loss low voltage": 0.0,
                "Transmission loss to low voltage": 0.0,
                "Production volume": 0.0,
            },
        )
        transf_loss += (
            dict_loss["Transformation loss low voltage"]
            * dict_loss["Production volume"]
        )
        distr_loss += (
            dict_loss["Transmission loss to low voltage"]
            * dict_loss["Production volume"]
        )
        cumul_prod += dict_loss["Production volume"]
    transf_loss /= cumul_prod
    distr_loss /= cumul_prod

    low = {"transf_loss": transf_loss, "distr_loss": distr_loss}

    return {"high": high, "medium": medium, "low": low}


class Electricity(BaseTransformation):
    """
    Class that modifies electricity markets in the database based on IAM output data.
    Inherits from `transformation.BaseTransformation`.

    :ivar database: wurst database, which is a list of dictionaries
    :vartype database: list
    :ivar iam_data: IAM data
    :vartype iam_data: xarray.DataArray
    :ivar model: name of the IAM model (e.g., "remind", "image")
    :vartype model: str
    :vartype pathway: str
    :ivar year: year of the pathway (e.g., 2030)
    :vartype year: int

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
    ) -> None:
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
        mapping = InventorySet(self.database)
        self.powerplant_map = mapping.generate_powerplant_map()
        # reverse dictionary of self.powerplant_map
        self.powerplant_map_rev = {}
        for k, v in self.powerplant_map.items():
            for pp in list(v):
                self.powerplant_map_rev[pp] = k

        self.powerplant_fuels_map = mapping.generate_powerplant_fuels_map()
        self.production_per_tech = self.get_production_per_tech_dict()
        losses = get_losses_per_country_dict()
        self.network_loss = {
            loc: get_production_weighted_losses(
                losses, self.geo.iam_to_ecoinvent_location(loc)
            )
            for loc in self.regions
        }
        self.system_model = system_model
        self.biosphere_dict = biosphere_flows_dictionary(self.version)

    @lru_cache
    def get_production_per_tech_dict(self) -> Dict[Tuple[str, str], float]:
        """
        Create a dictionary with tuples (technology, country) as keys
        and production volumes as values.
        :return: technology to production volume dictionary
        :rtype: dict
        """

        production_vols = {}

        for dataset_names in self.powerplant_map.values():
            for name in dataset_names:
                for dataset in ws.get_many(
                    self.database,
                    ws.equals("name", name),
                ):
                    for exc in ws.production(dataset):
                        # even if non-existent, we set a minimum value of 1e-9
                        # because if not, we risk dividing by zero!!!
                        production_vols[(dataset["name"], dataset["location"])] = max(
                            float(exc.get("production volume", 1e-9)), 1e-9
                        )

        return production_vols

    def check_for_production_volume(self, suppliers: List[dict]) -> List[dict]:
        # Remove suppliers that do not have a production volume
        return [
            supplier
            for supplier in suppliers
            if self.get_production_weighted_share(supplier, suppliers) != 0
        ]

    def get_production_weighted_share(
        self, supplier: dict, suppliers: List[dict]
    ) -> float:
        """
        Return the share of production of an electricity-producing dataset in a specific location,
        relative to the summed production of similar technologies in locations contained in
        the same IAM region.
        :param supplier: electricity-producing dataset
        :type supplier: wurst dataset
        :param suppliers: list of electricity-producing datasets
        :type suppliers: list of wurst datasets
        :return: share of production relative to the total population
        :rtype: float
        """

        # Fetch the production volume of the supplier
        loc_production = float(
            self.production_per_tech.get((supplier["name"], supplier["location"]), 0)
        )

        # Fetch the total production volume of similar technologies in other locations
        # contained within the IAM region.

        total_production = 0
        for loc in suppliers:
            total_production += float(
                self.production_per_tech.get((loc["name"], loc["location"]), 0)
            )

        # If a corresponding production volume is found.
        if total_production != 0:
            return loc_production / total_production
        else:
            # If not, we allocate an equal share of supply
            return 1 / len(suppliers)

    def create_new_markets_low_voltage(self) -> None:
        """
        Create low voltage market groups for electricity, by receiving medium voltage market groups as input
        and adding transformation and distribution losses. Transformation and distribution losses are taken from ei37.
        Contribution from solar power is added here as well, as most is delivered at low voltage,
        although CSP delivers at high voltage.
        Does not return anything. Modifies the database in place.
        """

        # Loop through the technologies
        technologies = [
            tech
            for tech in self.iam_data.electricity_markets.variables.values
            if "solar pv residential" in tech.lower()
        ]

        # Get the possible names of ecoinvent datasets
        ecoinvent_technologies = {
            technology: self.powerplant_map[technology] for technology in technologies
        }

        # Loop through IAM regions
        for region in self.regions:
            transf_loss = self.network_loss[region]["low"]["transf_loss"]
            distr_loss = self.network_loss[region]["low"]["distr_loss"]

            # Fetch ecoinvent regions contained in the IAM region
            ecoinvent_regions = self.geo.iam_to_ecoinvent_location(region)

            possible_locations = [
                [region],
                ecoinvent_regions,
                ["RER"],
                ["RoW"],
                ["CH"],
            ]

            tech_suppliers = defaultdict(list)

            for technology in ecoinvent_technologies:
                suppliers, counter = [], 0

                while len(suppliers) == 0:
                    suppliers = list(
                        get_suppliers_of_a_region(
                            database=self.database,
                            locations=possible_locations[counter],
                            names=ecoinvent_technologies[technology],
                            reference_prod="electricity",
                            unit="kilowatt hour",
                        )
                    )
                    counter += 1

                suppliers = self.check_for_production_volume(suppliers)
                for supplier in suppliers:
                    share = self.get_production_weighted_share(supplier, suppliers)
                    tech_suppliers[technology].append((supplier, share))

                # remove suppliers that have a supply share inferior to 0.1%
                tech_suppliers[technology] = [
                    supplier
                    for supplier in tech_suppliers[technology]
                    if supplier[1] > 0.001
                ]
                # rescale the shares so that they sum to 1
                total_share = sum(
                    [supplier[1] for supplier in tech_suppliers[technology]]
                )
                tech_suppliers[technology] = [
                    (supplier[0], supplier[1] / total_share)
                    for supplier in tech_suppliers[technology]
                ]

            # `period` is a period of time considered to create time-weighted average mix
            # when `period` == 0, this is a market mix for the year `self.year`
            # when `period` == 10, this is a market mix for the period `self.year` + 10
            # this is useful for systems that consume electricity
            # over a long period of time (e.g., buildings, BEVs, etc.)

            if self.system_model == "consequential":
                periods = [
                    0,
                ]
            else:
                periods = [0, 20, 40, 60]

            for period in periods:
                if self.system_model == "consequential":
                    electricity_mix = dict(
                        zip(
                            self.iam_data.electricity_markets.variables.values,
                            self.iam_data.electricity_markets.sel(
                                region=region, year=self.year
                            ).values,
                        )
                    )

                else:
                    # Create a time-weighted average mix
                    electricity_mix = dict(
                        zip(
                            self.iam_data.electricity_markets.variables.values,
                            self.iam_data.electricity_markets.sel(
                                region=region,
                            )
                            .interp(
                                year=np.arange(self.year, self.year + period + 1),
                                kwargs={"fill_value": "extrapolate"},
                            )
                            .mean(dim="year")
                            .values,
                        )
                    )

                # Create an empty dataset
                new_dataset = {
                    "location": region,
                    "name": "market group for electricity, low voltage",
                    "reference product": "electricity, low voltage",
                    "unit": "kilowatt hour",
                    "database": self.database[1]["database"],
                    "code": str(uuid.uuid4().hex),
                    "comment": f"Dataset created by `premise` from the IAM model {self.model.upper()}"
                    f" using the pathway {self.scenario} for the year {self.year}.",
                }

                # First, add the reference product exchange
                new_exchanges = [
                    {
                        "uncertainty type": 0,
                        "loc": 1,
                        "amount": 1,
                        "type": "production",
                        "production volume": 0,
                        "product": "electricity, low voltage",
                        "name": "market group for electricity, low voltage",
                        "unit": "kilowatt hour",
                        "location": region,
                    }
                ]

                if period != 0:
                    # this dataset is for a period of time
                    new_dataset["name"] += f", {period}-year period"
                    new_dataset["comment"] += (
                        f" Average electricity mix over a {period}"
                        f"-year period {self.year}-{self.year + period}."
                    )
                    new_exchanges[0]["name"] += f", {period}-year period"

                # Second, add an input of sulfur hexafluoride (SF6) emission to compensate the transformer's leakage
                # And an emission of a corresponding amount
                # Third, transmission line and SF6 supply and emission
                possible_suppliers = self.select_multiple_suppliers(
                    possible_names=("market for sulfur hexafluoride, liquid",),
                    dataset_location=region,
                )
                new_exchanges.extend(
                    [
                        {
                            "uncertainty type": 0,
                            "loc": supplier[1],
                            "amount": 2.99e-9 * share,
                            "type": "technosphere",
                            "production volume": 0,
                            "product": supplier[2],
                            "name": supplier[0],
                            "unit": supplier[-1],
                            "location": supplier[1],
                        }
                        for supplier, share in possible_suppliers.items()
                    ]
                )
                new_exchanges.append(
                    {
                        "uncertainty type": 0,
                        "loc": 2.99e-9,
                        "amount": 2.99e-9,
                        "type": "biosphere",
                        "input": (
                            "biosphere3",
                            self.biosphere_dict[
                                (
                                    "Sulfur hexafluoride",
                                    "air",
                                    "non-urban air or from high stacks",
                                    "kilogram",
                                )
                            ],
                        ),
                        "name": "Sulfur hexafluoride",
                        "unit": "kilogram",
                        "categories": ("air", "non-urban air or from high stacks"),
                    },
                )

                possible_suppliers = self.select_multiple_suppliers(
                    possible_names=(
                        "distribution network construction, electricity, low voltage",
                    ),
                    dataset_location=region,
                )
                new_exchanges.extend(
                    [
                        {
                            "uncertainty type": 0,
                            "loc": supplier[1],
                            "amount": 8.74e-8 * share,
                            "type": "technosphere",
                            "production volume": 0,
                            "product": supplier[2],
                            "name": supplier[0],
                            "unit": supplier[-1],
                            "location": supplier[1],
                        }
                        for supplier, share in possible_suppliers.items()
                    ]
                )

                # Fourth, add the contribution of solar power
                solar_amount = 0

                for technology in technologies:
                    # If the solar power technology contributes to the mix
                    if electricity_mix[technology] > 0:
                        # Contribution in supply
                        amount = electricity_mix[technology]
                        solar_amount += amount

                        for supplier, share in tech_suppliers[technology]:
                            new_exchanges.append(
                                {
                                    "uncertainty type": 0,
                                    "loc": (amount * share),
                                    "amount": (amount * share),
                                    "type": "technosphere",
                                    "production volume": 0,
                                    "product": supplier["reference product"],
                                    "name": supplier["name"],
                                    "unit": supplier["unit"],
                                    "location": supplier["location"],
                                }
                            )

                # Fifth, add:
                # * an input from the medium voltage market minus solar contribution, including distribution loss
                # * a self-consuming input for transformation loss

                new_exchanges.append(
                    {
                        "uncertainty type": 0,
                        "loc": 0,
                        "amount": (1 - solar_amount) * (1 + distr_loss),
                        "type": "technosphere",
                        "production volume": 0,
                        "product": "electricity, medium voltage",
                        "name": "market group for electricity, medium voltage"
                        if period == 0
                        else f"market group for electricity, medium voltage, {period}-year period",
                        "unit": "kilowatt hour",
                        "location": region,
                    }
                )

                new_exchanges.append(
                    {
                        "uncertainty type": 0,
                        "loc": 0,
                        "amount": transf_loss,
                        "type": "technosphere",
                        "production volume": 0,
                        "product": "electricity, low voltage",
                        "name": "market group for electricity, low voltage"
                        if period == 0
                        else f"market group for electricity, low voltage, {period}-year period",
                        "unit": "kilowatt hour",
                        "location": region,
                    }
                )

                new_dataset["exchanges"] = new_exchanges

                if "log parameters" not in new_dataset:
                    new_dataset["log parameters"] = {}

                new_dataset["log parameters"].update(
                    {
                        "distribution loss": distr_loss,
                        "transformation loss": transf_loss,
                        "renewable share": solar_amount / (1 + distr_loss),
                    }
                )
                self.database.append(new_dataset)

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

    def create_new_markets_medium_voltage(self) -> None:
        """
        Create medium voltage market groups for electricity, by receiving high voltage market groups as inputs
        and adding transformation and distribution losses.
        Contribution from solar power is added in low voltage market groups.
        Does not return anything. Modifies the database in place.
        """

        log_created_markets = []

        for region in self.regions:
            transf_loss = self.network_loss[region]["medium"]["transf_loss"]
            distr_loss = self.network_loss[region]["medium"]["distr_loss"]

            # `period` is a period of time considered to create time-weighted average mix
            # when `period` == 0, this is a market mix for the year `self.year`
            # when `period` == 10, this is a market mix for the period `self.year` + 10
            # this is useful for systems that consume electricity
            # over a long period of time (e.g., buildings, BEVs, etc.)

            if self.system_model == "consequential":
                periods = [
                    0,
                ]
            else:
                periods = [0, 20, 40, 60]

            for period in periods:
                # Create an empty dataset
                new_dataset = {
                    "location": region,
                    "name": "market group for electricity, medium voltage",
                    "reference product": "electricity, medium voltage",
                    "unit": "kilowatt hour",
                    "database": self.database[1]["database"],
                    "code": str(uuid.uuid4().hex),
                    "comment": f"Dataset created by `premise` from the IAM model {self.model.upper()}"
                    f" using the pathway {self.scenario} for the year {self.year}.",
                }

                # First, add the reference product exchange
                new_exchanges = [
                    {
                        "uncertainty type": 0,
                        "loc": 1,
                        "amount": 1,
                        "type": "production",
                        "production volume": 0,
                        "product": "electricity, medium voltage",
                        "name": "market group for electricity, medium voltage",
                        "unit": "kilowatt hour",
                        "location": region,
                    }
                ]

                if period != 0:
                    # this dataset is for a period of time
                    new_dataset["name"] += f", {period}-year period"
                    new_dataset["comment"] += (
                        f" Average electricity mix over a {period}"
                        f"-year period {self.year}-{self.year + period}."
                    )
                    new_exchanges[0]["name"] += f", {period}-year period"

                # Second, add:
                # * an input from the high voltage market, including transmission loss
                # * a self-consuming input for transformation loss

                new_exchanges.append(
                    {
                        "uncertainty type": 0,
                        "loc": 0,
                        "amount": 1 + distr_loss,
                        "type": "technosphere",
                        "production volume": 0,
                        "product": "electricity, high voltage",
                        "name": "market group for electricity, high voltage"
                        if period == 0
                        else f"market group for electricity, high voltage, {period}-year period",
                        "unit": "kilowatt hour",
                        "location": region,
                    }
                )

                new_exchanges.append(
                    {
                        "uncertainty type": 0,
                        "loc": 0,
                        "amount": transf_loss,
                        "type": "technosphere",
                        "production volume": 0,
                        "product": "electricity, medium voltage",
                        "name": "market group for electricity, medium voltage"
                        if period == 0
                        else f"market group for electricity, medium voltage, {period}-year period",
                        "unit": "kilowatt hour",
                        "location": region,
                    }
                )

                # Third, add an input to of sulfur hexafluoride emission to compensate the transformer's leakage
                # And an emission of a corresponding amount

                possible_suppliers = self.select_multiple_suppliers(
                    possible_names=("market for sulfur hexafluoride, liquid",),
                    dataset_location=region,
                )
                new_exchanges.extend(
                    [
                        {
                            "uncertainty type": 0,
                            "loc": supplier[1],
                            "amount": 5.4e-8 * share,
                            "type": "technosphere",
                            "production volume": 0,
                            "product": supplier[2],
                            "name": supplier[0],
                            "unit": supplier[-1],
                            "location": supplier[1],
                        }
                        for supplier, share in possible_suppliers.items()
                    ]
                )
                new_exchanges.append(
                    {
                        "uncertainty type": 0,
                        "loc": 5.4e-8,
                        "amount": 5.4e-8,
                        "type": "biosphere",
                        "input": (
                            "biosphere3",
                            self.biosphere_dict[
                                (
                                    "Sulfur hexafluoride",
                                    "air",
                                    "non-urban air or from high stacks",
                                    "kilogram",
                                )
                            ],
                        ),
                        "name": "Sulfur hexafluoride",
                        "unit": "kilogram",
                        "categories": ("air", "non-urban air or from high stacks"),
                    },
                )

                # Fourth, transmission line
                possible_suppliers = self.select_multiple_suppliers(
                    possible_names=(
                        "transmission network construction, electricity, medium voltage",
                    ),
                    dataset_location=region,
                )
                new_exchanges.extend(
                    [
                        {
                            "uncertainty type": 0,
                            "loc": supplier[1],
                            "amount": 1.8628e-8 * share,
                            "type": "technosphere",
                            "production volume": 0,
                            "product": supplier[2],
                            "name": supplier[0],
                            "unit": supplier[-1],
                            "location": supplier[1],
                        }
                        for supplier, share in possible_suppliers.items()
                    ]
                )

                new_dataset["exchanges"] = new_exchanges

                if "log parameters" not in new_dataset:
                    new_dataset["log parameters"] = {}

                new_dataset["log parameters"].update(
                    {
                        "distribution loss": distr_loss,
                        "transformation loss": transf_loss,
                        "renewable share": 0.0,
                    }
                )

                self.database.append(new_dataset)

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

    def create_new_markets_high_voltage(self) -> None:
        """
        Create high voltage market groups for electricity, based on electricity mixes given by the IAM scenario.
        Contribution from solar power is added in low voltage market groups.
        Does not return anything. Modifies the database in place.
        """

        log_created_markets = []

        # Loop through the technologies
        technologies = [
            tech
            for tech in self.iam_data.electricity_markets.variables.values
            if "solar pv residential" not in tech.lower()
        ]

        # Get the possible names of ecoinvent datasets
        ecoinvent_technologies = {
            technology: self.powerplant_map[technology] for technology in technologies
        }

        for region in self.regions:
            # Fetch ecoinvent regions contained in the IAM region
            ecoinvent_regions = self.geo.iam_to_ecoinvent_location(region)
            # Second, add transformation loss
            transf_loss = self.network_loss[region]["high"]["transf_loss"]

            # Fetch electricity-producing technologies contained in the IAM region
            # if they cannot be found for the ecoinvent locations concerned
            # we widen the scope to EU-based datasets, and RoW, and finally Switzerland

            possible_locations = [
                [region],
                ecoinvent_regions,
                ["RER"],
                ["RoW"],
                ["CH"],
            ]

            tech_suppliers = defaultdict(list)

            for technology in ecoinvent_technologies:
                suppliers, counter = [], 0

                try:
                    while len(suppliers) == 0:
                        suppliers = list(
                            get_suppliers_of_a_region(
                                database=self.database,
                                locations=possible_locations[counter],
                                names=ecoinvent_technologies[technology],
                                reference_prod="electricity",
                                unit="kilowatt hour",
                            )
                        )
                        counter += 1

                    suppliers = self.check_for_production_volume(suppliers)

                    for supplier in suppliers:
                        share = self.get_production_weighted_share(supplier, suppliers)
                        tech_suppliers[technology].append((supplier, share))

                    # remove suppliers that have a supply share inferior to 0.1%
                    tech_suppliers[technology] = [
                        supplier
                        for supplier in tech_suppliers[technology]
                        if supplier[1] > 0.001
                    ]
                    # rescale the shares so that they sum to 1
                    total_share = sum(
                        [supplier[1] for supplier in tech_suppliers[technology]]
                    )
                    tech_suppliers[technology] = [
                        (supplier[0], supplier[1] / total_share)
                        for supplier in tech_suppliers[technology]
                    ]

                except IndexError:
                    if self.system_model == "consequential":
                        continue
                    else:
                        raise

            if self.system_model == "consequential":
                periods = [
                    0,
                ]
            else:
                periods = [0, 20, 40, 60]

            for period in periods:
                if self.system_model == "consequential":
                    electricity_mix = dict(
                        zip(
                            self.iam_data.electricity_markets.variables.values,
                            self.iam_data.electricity_markets.sel(
                                region=region, year=self.year
                            ).values,
                        )
                    )

                else:
                    electricity_mix = dict(
                        zip(
                            self.iam_data.electricity_markets.variables.values,
                            self.iam_data.electricity_markets.sel(
                                region=region,
                            )
                            .interp(
                                year=np.arange(self.year, self.year + period + 1),
                                kwargs={"fill_value": "extrapolate"},
                            )
                            .mean(dim="year")
                            .values,
                        )
                    )

                new_dataset = {
                    "location": region,
                    "name": "market group for electricity, high voltage",
                    "reference product": "electricity, high voltage",
                    "unit": "kilowatt hour",
                    "database": self.database[1]["database"],
                    "code": str(uuid.uuid4().hex),
                    "comment": f"Dataset created by `premise` from the IAM model {self.model.upper()}"
                    f" using the pathway {self.scenario} for the year {self.year}.",
                }

                # First, add the reference product exchange
                new_exchanges = [
                    {
                        "uncertainty type": 0,
                        "loc": 1,
                        "amount": 1,
                        "type": "production",
                        "production volume": 0,
                        "product": "electricity, high voltage",
                        "name": "market group for electricity, high voltage",
                        "unit": "kilowatt hour",
                        "location": region,
                    }
                ]

                new_exchanges.append(
                    {
                        "uncertainty type": 0,
                        "loc": 1,
                        "amount": transf_loss,
                        "type": "technosphere",
                        "production volume": 0,
                        "product": "electricity, high voltage",
                        "name": "market group for electricity, high voltage",
                        "unit": "kilowatt hour",
                        "location": region,
                    }
                )

                if period != 0:
                    # this dataset is for a period of time
                    new_dataset["name"] += f", {period}-year period"
                    new_dataset["comment"] += (
                        f" Average electricity mix over a {period}"
                        f"-year period {self.year}-{self.year + period}."
                    )
                    new_exchanges[0]["name"] += f", {period}-year period"
                    new_exchanges[-1]["name"] += f", {period}-year period"

                # Fetch residential solar PV contribution in the mix, to subtract it
                # as solar energy is an input of low-voltage markets
                solar_amount = 0
                for tech in electricity_mix:
                    if "solar pv residential" in tech.lower():
                        solar_amount += electricity_mix[tech]

                # calculate the share of renewable energy in the mix
                renewable_share = 0
                renewable_techs = [
                    "solar",
                    "wind",
                    "geothermal",
                    "hydro",
                    "biomass",
                    "biogas",
                    "wave",
                ]
                for tech in electricity_mix:
                    if any(x in tech.lower() for x in renewable_techs):
                        renewable_share += electricity_mix[tech]

                for technology in technologies:
                    # If the given technology contributes to the mix
                    if electricity_mix[technology] > 0:
                        # Contribution in supply
                        amount = electricity_mix[technology]

                        for supplier, share in tech_suppliers[technology]:
                            new_exchanges.append(
                                {
                                    "uncertainty type": 0,
                                    "loc": (amount * share) / (1 - solar_amount),
                                    "amount": (amount * share) / (1 - solar_amount),
                                    "type": "technosphere",
                                    "production volume": 0,
                                    "product": supplier["reference product"],
                                    "name": supplier["name"],
                                    "unit": supplier["unit"],
                                    "location": supplier["location"],
                                }
                            )

                            log_created_markets.append(
                                [
                                    f"high voltage, {self.scenario}, {self.year}"
                                    if period == 0
                                    else f"high voltage, {self.scenario}, {self.year}, {period}-year period",
                                    technology,
                                    region,
                                    transf_loss,
                                    0.0,
                                    supplier["name"],
                                    supplier["location"],
                                    share,
                                    (amount * share) / (1 - solar_amount),
                                ]
                            )

                new_dataset["exchanges"] = new_exchanges

                if "log parameters" not in new_dataset:
                    new_dataset["log parameters"] = {}

                new_dataset["log parameters"].update(
                    {
                        "distribution loss": 0.0,
                        "transformation loss": transf_loss,
                        "renewable share": (renewable_share - solar_amount)
                        / sum(electricity_mix.values()),
                    }
                )

                self.database.append(new_dataset)

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

    def update_efficiency_of_solar_pv(self) -> None:
        """
        Update the efficiency of solar PV modules.
        We look at how many square meters are needed per kilowatt of installed capacity
        to obtain the current efficiency.
        Then we update the surface needed according to the projected efficiency.
        :return:
        """

        print("Update efficiency of solar PV panels.")

        # TODO: check if IAM data provides efficiencies for PV panels and use them instead

        # efficiency of modules in the future
        module_eff = get_efficiency_solar_photovoltaics()

        datasets = ws.get_many(
            self.database,
            *[
                ws.contains("name", "photovoltaic"),
                ws.either(
                    ws.contains("name", "installation"),
                    ws.contains("name", "construction"),
                ),
                ws.doesnt_contain_any("name", ["market", "factory", "module"]),
                ws.equals("unit", "unit"),
            ],
        )

        for dataset in datasets:
            power = float(re.findall(r"[-+]?\d*\.\d+|\d+", dataset["name"])[0])

            if "mwp" in dataset["name"].lower():
                power *= 1000

            for exc in ws.technosphere(
                dataset,
                *[
                    ws.contains("name", "photovoltaic"),
                    ws.equals("unit", "square meter"),
                ],
            ):
                surface = float(exc["amount"])
                max_power = surface  # in kW, since we assume a constant 1,000W/m^2
                current_eff = power / max_power

                possible_techs = [
                    "micro-Si",
                    "single-Si",
                    "multi-Si",
                    "CIGS",
                    "CIS",
                    "CdTe",
                ]
                pv_tech = [
                    i for i in possible_techs if i.lower() in exc["name"].lower()
                ]

                if len(pv_tech) > 0:
                    pv_tech = pv_tech[0]

                    new_eff = (
                        module_eff.sel(technology=pv_tech)
                        .interp(year=self.year, kwargs={"fill_value": "extrapolate"})
                        .values
                    )

                    # in case self.year <10 or >2050
                    new_eff = np.clip(new_eff, 0.1, 0.27)

                    # We only update the efficiency if it is higher than the current one.
                    if new_eff > current_eff:
                        exc["amount"] *= float(current_eff / new_eff)

                        dataset["comment"] = (
                            f"`premise` has changed the efficiency "
                            f"of this photovoltaic installation "
                            f"from {int(current_eff * 100)} pct. to {int(new_eff * 100)} pt."
                        )

                        if "log parameters" not in dataset:
                            dataset["log parameters"] = {}

                        dataset["log parameters"].update(
                            {"old efficiency": current_eff}
                        )
                        dataset["log parameters"].update({"new efficiency": new_eff})

                        # add to log
                        self.write_log(dataset=dataset, status="updated")

    def update_ng_production_ds(self) -> None:
        """
        Relink updated datasets for natural gas extraction from
        http://www.esu-services.ch/fileadmin/download/publicLCI/meili-2021-LCI%20for%20the%20oil%20and%20gas%20extraction.pdf
        to high pressure natural gas markets.
        """

        print("Update natural gas extraction datasets.")

        countries = ["NL", "DE", "FR", "RER", "IT", "CH"]

        for dataset in self.database:
            amount = {}
            to_remove = []
            for exc in dataset["exchanges"]:
                if (
                    exc["name"] == "market for natural gas, high pressure"
                    and exc["location"] in countries
                    and exc["type"] == "technosphere"
                ):
                    if exc["location"] in amount:
                        amount[exc["location"]] += exc["amount"]
                    else:
                        amount[exc["location"]] = exc["amount"]
                    to_remove.append(
                        (exc["name"], exc["product"], exc["location"], exc["type"])
                    )

            if amount:
                dataset["exchanges"] = [
                    e
                    for e in dataset["exchanges"]
                    if (e["name"], e.get("product"), e.get("location"), e["type"])
                    not in to_remove
                ]

                for loc in amount:
                    dataset["exchanges"].append(
                        {
                            "name": "natural gas, high pressure, at consumer",
                            "product": "natural gas, high pressure, at consumer",
                            "location": loc,
                            "unit": "cubic meter",
                            "amount": amount[loc],
                            "type": "technosphere",
                        }
                    )

        countries = ["DE", "DZ", "GB", "NG", "NL", "NO", "RU", "US"]

        names = ["natural gas production", "petroleum and gas production"]

        for dataset in self.database:
            amount = {}
            to_remove = []
            for exc in dataset["exchanges"]:
                if (
                    any(i in exc["name"] for i in names)
                    and "natural gas, high pressure"
                    and exc["location"] in countries
                    and exc["type"] == "technosphere"
                ):
                    if exc["location"] in amount:
                        amount[exc["location"]] += exc["amount"]
                    else:
                        amount[exc["location"]] = exc["amount"]
                    to_remove.append(
                        (exc["name"], exc["product"], exc["location"], exc["type"])
                    )

            if amount:
                dataset["exchanges"] = [
                    e
                    for e in dataset["exchanges"]
                    if (e["name"], e.get("product"), e.get("location"), e["type"])
                    not in to_remove
                ]

                for loc in amount:
                    dataset["exchanges"].append(
                        {
                            "name": "natural gas, at production",
                            "product": "natural gas, high pressure",
                            "location": loc,
                            "unit": "cubic meter",
                            "amount": amount[loc],
                            "type": "technosphere",
                        }
                    )

    def create_biomass_markets(self) -> None:
        print("Create biomass markets.")

        with open(IAM_BIOMASS_VARS, "r", encoding="utf-8") as stream:
            biomass_map = yaml.safe_load(stream)

        # create region-specific "Supply of forest residue" datasets
        forest_residues_ds = self.fetch_proxies(
            name=biomass_map["biomass - residual"]["ecoinvent_aliases"]["name"][0],
            ref_prod=biomass_map["biomass - residual"]["ecoinvent_aliases"][
                "reference product"
            ][0],
            production_variable="biomass - residual",
        )

        # add them to the database
        self.database.extend(forest_residues_ds.values())

        # add log
        for dataset in list(forest_residues_ds.values()):
            self.write_log(dataset=dataset)
            # add it to list of created datasets
            self.modified_datasets[(self.model, self.scenario, self.year)][
                "created"
            ].append(
                (
                    dataset["name"],
                    dataset["reference product"],
                    dataset["location"],
                    dataset["unit"],
                )
            )

        for region in self.regions:
            dataset = {
                "name": "market for biomass, used as fuel",
                "reference product": "biomass, used as fuel",
                "location": region,
                "comment": f"Biomass market, created by `premise`, "
                f"to align with projections for the region {region} in {self.year}. "
                "Calculated for an average energy input (LHV) of 19 MJ/kg, dry basis. "
                "Sum of inputs can be superior to 1, as "
                "inputs of wood chips, wet-basis, have been multiplied by a factor 2.5, "
                "to reach a LHV of 19 MJ (they have a LHV of 7.6 MJ, wet basis).",
                "unit": "kilogram",
                "database": eidb_label(
                    self.model,
                    self.scenario,
                    self.year,
                    self.version,
                    self.system_model,
                ),
                "code": str(uuid.uuid4().hex),
                "exchanges": [
                    {
                        "name": "market for biomass, used as fuel",
                        "product": "biomass, used as fuel",
                        "amount": 1,
                        "unit": "kilogram",
                        "location": region,
                        "uncertainty type": 0,
                        "type": "production",
                    }
                ],
            }

            for biomass_type, biomass_act in biomass_map.items():
                total_prod_vol = np.clip(
                    (
                        self.iam_data.production_volumes.sel(
                            variables=list(biomass_map.keys()), region=region
                        )
                        .interp(year=self.year)
                        .sum(dim="variables")
                    ),
                    1e-6,
                    None,
                )

                share = np.clip(
                    (
                        self.iam_data.production_volumes.sel(
                            variables=biomass_type, region=region
                        )
                        .interp(year=self.year)
                        .sum()
                        / total_prod_vol
                    ).values.item(0),
                    0,
                    1,
                )

                if not share:
                    if biomass_type == "biomass - residual":
                        share = 1
                    else:
                        share = 0

                if share > 0:
                    ecoinvent_regions = self.geo.iam_to_ecoinvent_location(
                        dataset["location"]
                    )
                    possible_locations = [
                        dataset["location"],
                        *ecoinvent_regions,
                        "RER",
                        "Europe without Switzerland",
                        "RoW",
                        "GLO",
                    ]
                    possible_names = biomass_act["ecoinvent_aliases"]["name"]
                    possible_products = biomass_act["ecoinvent_aliases"][
                        "reference product"
                    ]

                    suppliers, counter = [], 0

                    while not suppliers:
                        suppliers = list(
                            ws.get_many(
                                self.database,
                                ws.either(
                                    *[
                                        ws.contains("name", sup)
                                        for sup in possible_names
                                    ]
                                ),
                                ws.equals("location", possible_locations[counter]),
                                ws.either(
                                    *[
                                        ws.contains("reference product", prod)
                                        for prod in possible_products
                                    ]
                                ),
                                ws.equals("unit", "kilogram"),
                                ws.doesnt_contain_any(
                                    "name", ["willow", "post-consumer"]
                                ),
                            )
                        )
                        counter += 1

                    suppliers = get_shares_from_production_volume(suppliers)

                    for supplier, supply_share in suppliers.items():
                        multiplication_factor = 1.0
                        amount = supply_share * share * multiplication_factor
                        dataset["exchanges"].append(
                            {
                                "type": "technosphere",
                                "product": supplier[2],
                                "name": supplier[0],
                                "unit": supplier[-1],
                                "location": supplier[1],
                                "amount": amount,
                                "uncertainty type": 0,
                            }
                        )

                if "log parameters" not in dataset:
                    dataset["log parameters"] = {}

                dataset["log parameters"].update(
                    {
                        "biomass share": share,
                    }
                )

            self.database.append(dataset)

            # add log
            self.write_log(dataset=dataset)

            # add it to list of created datasets
            self.modified_datasets[(self.model, self.scenario, self.year)][
                "created"
            ].append(
                (
                    dataset["name"],
                    dataset["reference product"],
                    dataset["location"],
                    dataset["unit"],
                )
            )

        # replace biomass inputs
        print("Replace biomass inputs.")
        for dataset in ws.get_many(
            self.database,
            ws.either(
                *[ws.equals("unit", unit) for unit in ["kilowatt hour", "megajoule"]]
            ),
            ws.either(
                *[
                    ws.contains("name", name)
                    for name in ["electricity", "heat", "power"]
                ]
            ),
        ):
            for exc in ws.technosphere(
                dataset,
                ws.contains("name", "market for wood chips"),
                ws.equals("unit", "kilogram"),
            ):
                exc["name"] = "market for biomass, used as fuel"
                exc["product"] = "biomass, used as fuel"
                exc["location"] = self.ecoinvent_to_iam_loc[dataset["location"]]

        mapping = InventorySet(self.database)
        self.powerplant_fuels_map = mapping.generate_powerplant_fuels_map()

    def create_region_specific_power_plants(self):
        """
        Some power plant inventories are not native to ecoinvent
        but imported. However, they are defined for a specific location
        (mostly European), but are used in many electricity markets
        (non-European). Hence, we create region-specific versions of these datasets,
        to align inputs providers with the geographical scope of the region.

        """

        print("Create region-specific power plants.")
        all_plants = []

        techs = [
            "Biomass CHP CCS",
            "Biomass ST",
            "Biomass IGCC CCS",
            "Biomass IGCC",
            "Coal IGCC",
            "Coal PC CCS",
            "Coal CHP CCS",
            "Coal IGCC CCS",
            "Gas CHP CCS",
            "Gas CC CCS",
            "Oil CC CCS",
        ]

        list_datasets_to_duplicate = [
            dataset["name"]
            for dataset in self.database
            if dataset["name"]
            in [y for k, v in self.powerplant_map.items() for y in v if k in techs]
        ]

        list_datasets_to_duplicate.extend(
            [
                "Wood chips, burned in power plant",
                "Natural gas, in ATR ",
                "Hard coal, burned",
                "Lignite, burned",
                "CO2 storage/",
                "CO2 capture/",
            ]
        )

        for dataset in ws.get_many(
            self.database,
            ws.either(
                *[ws.contains("name", name) for name in list_datasets_to_duplicate]
            ),
        ):
            new_plants = self.fetch_proxies(
                name=dataset["name"],
                ref_prod=dataset["reference product"],
                production_variable=self.powerplant_map_rev.get(dataset["name"]),
            )

            # we need to adjust the need to CO2 capture and storage
            # based on the electricity provider in the dataset
            # hence, we want to know how much CO2 is released
            # by each provider, and capture 90% of the amount

            if "CHP CCS" in self.powerplant_map_rev.get(dataset["name"], ""):
                for plant in new_plants.values():
                    co2_amount = 0

                    providers = [
                        e
                        for e in plant["exchanges"]
                        if e["type"] == "technosphere" and e["unit"] == "kilowatt hour"
                    ]

                    for provider in providers:
                        provider_ds = ws.get_one(
                            self.database,
                            ws.equals("name", provider["name"]),
                            ws.equals("location", provider["location"]),
                            ws.equals("reference product", provider["product"]),
                            ws.equals("unit", provider["unit"]),
                        )
                        co2_amount += sum(
                            f["amount"] * provider["amount"]
                            for f in ws.biosphere(
                                provider_ds,
                                ws.contains("name", "Carbon dioxide"),
                            )
                        )

                    for exc in plant["exchanges"]:
                        if (
                            exc["type"] == "technosphere"
                            and exc["unit"] == "kilogram"
                            and exc["name"].startswith("CO2 capture")
                        ):
                            exc["amount"] = co2_amount * 0.9

                        if (
                            exc["type"] == "biosphere"
                            and exc["unit"] == "kilogram"
                            and exc["name"].startswith("Carbon dioxide")
                        ):
                            exc["amount"] = co2_amount * 0.9

            all_plants.extend(new_plants.values())

        self.database.extend(all_plants)

        for dataset in all_plants:
            self.write_log(dataset=dataset)
            # add it to list of created datasets
            self.modified_datasets[(self.model, self.scenario, self.year)][
                "created"
            ].append(
                (
                    dataset["name"],
                    dataset["reference product"],
                    dataset["location"],
                    dataset["unit"],
                )
            )

    def update_electricity_efficiency(self) -> None:
        """
        This method modifies each ecoinvent coal, gas,
        oil and biomass dataset using data from the IAM scenario.
        Return a wurst database with modified datasets.

        :return: a wurst database, with rescaled electricity-producing datasets.
        :rtype: list
        """

        print("Adjust efficiency of power plants...")

        mapping = InventorySet(self.database)
        self.fuel_map = mapping.generate_fuel_map()
        # reverse the fuel map to get a mapping from ecoinvent to premise
        self.fuel_map_reverse: Dict = {}

        for key, value in self.fuel_map.items():
            for v in list(value):
                self.fuel_map_reverse[v] = key

        eff_labels = self.iam_data.electricity_efficiencies.variables.values
        all_techs = self.iam_data.electricity_markets.variables.values

        technologies_map = self.get_iam_mapping(
            activity_map=self.powerplant_map,
            fuels_map=self.powerplant_fuels_map,
            technologies=list(set(eff_labels).intersection(all_techs)),
        )

        # to store changes in efficiency
        eff_change_log = []

        for technology in technologies_map:
            dict_technology = technologies_map[technology]
            print("Rescale inventories and emissions for", technology)

            for dataset in ws.get_many(
                self.database,
                ws.equals("unit", "kilowatt hour"),
                ws.either(
                    *[
                        ws.contains("name", n)
                        for n in dict_technology["technology filters"]
                    ]
                ),
            ):
                # Find current efficiency
                ei_eff = dict_technology["current_eff_func"](
                    dataset, dict_technology["fuel filters"], 3.6
                )

                # Find relative efficiency change indicated by the IAM
                scaling_factor = 1 / self.find_iam_efficiency_change(
                    data=self.iam_data.electricity_efficiencies,
                    variable=technology,
                    location=self.geo.ecoinvent_to_iam_location(dataset["location"]),
                )

                new_efficiency = float(np.clip(ei_eff * 1 / scaling_factor, 0, 1))

                if "log parameters" not in dataset:
                    dataset["log parameters"] = {}

                dataset["log parameters"].update(
                    {
                        "old efficiency": ei_eff,
                        "new efficiency": new_efficiency,
                    }
                )

                self.update_ecoinvent_efficiency_parameter(
                    dataset, ei_eff, new_efficiency
                )

                # Rescale all the technosphere exchanges
                # according to the change in efficiency between `year`
                # and 2020 from the IAM efficiency values
                wurst.change_exchanges_by_constant_factor(
                    dataset,
                    scaling_factor,
                )

                self.write_log(dataset=dataset, status="updated")

    def update_electricity_markets(self) -> None:
        """
        Delete electricity markets. Create high, medium and low voltage market groups for electricity.
        Link electricity-consuming datasets to newly created market groups for electricity.
        Return a wurst database with modified datasets.

        :return: a wurst database with new market groups for electricity
        :rtype: list
        """

        list_to_empty = [
            "market group for electricity",
            "market for electricity",
            "electricity, high voltage, import",
            "electricity, high voltage, production mix",
        ]

        # we want to preserve some electricity-related datasets
        list_to_preserve = [
            "cobalt industry",
            "aluminium industry",
            "coal mining",
            "label-certified",
            "renewable energy products",
            "for reuse in municipal waste incineration",
            "Swiss Federal Railways",
        ]

        # We first need to empty 'market for electricity' and 'market group for electricity' datasets
        print("Empty old electricity datasets")

        datasets_to_empty = ws.get_many(
            self.database,
            ws.either(*[ws.contains("name", n) for n in list_to_empty]),
            ws.equals("unit", "kilowatt hour"),
            ws.doesnt_contain_any("name", list_to_preserve),
        )

        for dataset in datasets_to_empty:
            # add tag
            dataset["has_downstream_consumer"] = False
            dataset["exchanges"] = [
                e for e in dataset["exchanges"] if e["type"] == "production"
            ]

            if "high voltage" in dataset["name"]:
                voltage = "high voltage"
            elif "medium voltage" in dataset["name"]:
                voltage = "medium voltage"
            else:
                voltage = "low voltage"

            dataset["exchanges"].append(
                {
                    "name": f"market group for electricity, {voltage}",
                    "product": f"electricity, {voltage}",
                    "amount": 1,
                    "uncertainty type": 0,
                    "location": self.ecoinvent_to_iam_loc[dataset["location"]],
                    "type": "technosphere",
                    "unit": "kilowatt hour",
                }
            )

            self.write_log(dataset=dataset, status="updated")

            # list `market group for electricity` as "emptied"
            self.modified_datasets[(self.model, self.scenario, self.year)][
                "emptied"
            ].append(
                (
                    dataset["name"],
                    dataset["reference product"],
                    dataset["location"],
                    dataset["unit"],
                )
            )

            # add new regional datasets to cache
            self.add_new_entry_to_cache(
                location=dataset["location"],
                exchange=dataset,
                allocated=[
                    {
                        "name": f"market group for electricity, {voltage}",
                        "reference product": f"electricity, {voltage}",
                        "unit": "kilowatt hour",
                        "location": self.ecoinvent_to_iam_loc[dataset["location"]],
                    },
                ],
                shares=[
                    1.0,
                ],
            )

        # We then need to create high voltage IAM electricity markets
        print("Create high voltage markets.")
        self.create_new_markets_high_voltage()
        print("Create medium voltage markets.")
        self.create_new_markets_medium_voltage()
        print("Create low voltage markets.")
        self.create_new_markets_low_voltage()

        print("Done!")

    def write_log(self, dataset, status="created"):
        """
        Write log file.
        """

        logger.info(
            f"{status}|{self.model}|{self.scenario}|{self.year}|"
            f"{dataset['name']}|{dataset['location']}|"
            f"{dataset.get('log parameters', {}).get('old efficiency', '')}|"
            f"{dataset.get('log parameters', {}).get('new efficiency', '')}|"
            f"{dataset.get('log parameters', {}).get('biomass share', '')}|"
            f"{dataset.get('log parameters', {}).get('transformation loss', '')}|"
            f"{dataset.get('log parameters', {}).get('distribution loss', '')}|"
            f"{dataset.get('log parameters', {}).get('renewable share', '')}"
        )
