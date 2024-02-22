"""
Integrates projections regarding fuel production and supply.
"""

import copy
from functools import lru_cache
from typing import Union

import xarray as xr
import yaml
from numpy import ndarray

from .filesystem_constants import DATA_DIR, VARIABLES_DIR
from .inventory_imports import get_biosphere_code
from .logger import create_logger
from .transformation import (
    Any,
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
from .utils import get_crops_properties, rescale_exchanges

logger = create_logger("fuel")

REGION_CLIMATE_MAP = VARIABLES_DIR / "iam_region_to_climate.yaml"
FUEL_LABELS = DATA_DIR / "fuels" / "fuel_labels.csv"
SUPPLY_CHAIN_SCENARIOS = DATA_DIR / "fuels" / "supply_chain_scenarios.yml"
HEAT_SOURCES = DATA_DIR / "fuels" / "heat_sources_map.yml"
HYDROGEN_SOURCES = DATA_DIR / "fuels" / "hydrogen_efficiency_parameters.yml"
HYDROGEN_SUPPLY_LOSSES = DATA_DIR / "fuels" / "hydrogen_supply_losses.yml"
METHANE_SOURCES = DATA_DIR / "fuels" / "methane_activities.yml"
LIQUID_FUEL_SOURCES = DATA_DIR / "fuels" / "liquid_fuel_activities.yml"
FUEL_MARKETS = DATA_DIR / "fuels" / "fuel_markets.yml"
BIOFUEL_SOURCES = DATA_DIR / "fuels" / "biofuels_activities.yml"
FUEL_GROUPS = DATA_DIR / "fuels" / "fuel_groups.yaml"


def fetch_mapping(filepath: str) -> dict:
    """Returns a dictionary from a YML file"""

    with open(filepath, encoding="utf-8") as stream:
        mapping = yaml.safe_load(stream)
    return mapping


@lru_cache()
def get_compression_effort(
    inlet_pressure: int, outlet_pressure: int, flow_rate: int
) -> float:
    """
    Calculate the required electricity consumption from the compressor given
    an inlet and outlet pressure and a flow rate for hydrogen.

    :param inlet_pressure: the input pressure in bar
    :param outlet_pressure: the output pressure in bar
    :param flow_rate: the flow rate of hydrogen in kg/day
    :return: the required electricity consumption in kWh

    """
    # Constants
    COMPRESSIBILITY_FACTOR = 1.03198
    NUM_COMPRESSOR_STAGES = 2
    INLET_TEMPERATURE = 310.95  # K
    RATIO_SPECIFIC_HEATS = 1.4
    MOLECULAR_MASS_H2 = 2.15  # g/mol
    COMPRESSOR_EFFICIENCY = 0.75

    # Intermediate calculations
    mass_flow_rate = flow_rate / (24 * 3600)  # convert to kg/s
    specific_gas_constant = 8.314  # J/(mol*K)
    part_1 = (
        mass_flow_rate
        * (COMPRESSIBILITY_FACTOR * INLET_TEMPERATURE * specific_gas_constant)
        / (MOLECULAR_MASS_H2 * COMPRESSOR_EFFICIENCY)
        * (NUM_COMPRESSOR_STAGES * RATIO_SPECIFIC_HEATS / (RATIO_SPECIFIC_HEATS - 1))
    )
    part_2 = (
        (outlet_pressure / inlet_pressure)
        ** ((RATIO_SPECIFIC_HEATS - 1) / (NUM_COMPRESSOR_STAGES * RATIO_SPECIFIC_HEATS))
    ) - 1
    shaft_power = part_1 * part_2

    # Convert to kWh
    electricity_consumption = shaft_power * 24 / 1000

    return electricity_consumption


@lru_cache()
def get_pre_cooling_energy(
    ambient_temperature: float, capacity_utilization: float
) -> float:
    """
    Calculate the required electricity consumption to pre-cool the hydrogen
    before tank filling.

    :param ambient_temperature: the ambient temperature in degrees Celsius
    :param capacity_utilization: the capacity utilization of the pre-cooling system
    :return: the required electricity consumption in kWh

    """
    # Constants
    COEFFICIENT_1 = 0.3 / 1.6
    COEFFICIENT_2 = -0.018
    COEFFICIENT_3 = 25
    COEFFICIENT_4 = -21

    # Convert temperature to Kelvin
    temperature_K = ambient_temperature + 273.15

    # Calculate pre-cooling energy
    energy_pre_cooling = (
        COEFFICIENT_1 * np.exp(COEFFICIENT_2 * ambient_temperature)
        + (COEFFICIENT_3 * np.log(temperature_K) + COEFFICIENT_4) / capacity_utilization
    )

    return energy_pre_cooling


@lru_cache()
def adjust_electrolysis_electricity_requirement(year: int) -> ndarray:
    """

    Calculate the adjusted electricity requirement for hydrogen electrolysis
    based on the given year.

    The electricity requirement decreases linearly from 58 kWh/kg H2 in 2010
    to 48 kWh/kg H2 in 2050, according to a literature review conducted by
    the Paul Scherrer Institute:

    Bauer (ed.), C., Desai, H., Heck, T., Sacchi, R., Schneider, S., Terlouw,
    T., Treyer, K., Zhang, X. Electricity storage and hydrogen – technologies,
    costs and impacts on climate change.
    Auftraggeberin: Bundesamt für Energie BFE, 3003 Bern.


    :param year: the year for which to calculate the adjusted electricity requirement
    :return: the adjusted electricity requirement in kWh/kg H2

    """
    # Constants
    MIN_ELECTRICITY_REQUIREMENT = 48
    MAX_ELECTRICITY_REQUIREMENT = 60  # no maximum

    # Calculate adjusted electricity requirement
    electricity_requirement = -0.3538 * (year - 2010) + 58.589

    # Clip to minimum and maximum values
    adjusted_requirement = np.clip(
        electricity_requirement,
        MIN_ELECTRICITY_REQUIREMENT,
        MAX_ELECTRICITY_REQUIREMENT,
    )

    return adjusted_requirement


def is_fuel_production(name):
    return any(
        i.lower() in name.lower()
        for i in ["ethanol production", "biodiesel production"]
    )


def update_co2_emissions(
    dataset: dict, amount_non_fossil_co2: float, biosphere_flows: dict
) -> dict:
    """Update fossil and non-fossil CO2 emissions of the dataset."""
    # Test for the presence of a fossil CO2 flow
    if not any(
        exc for exc in dataset["exchanges"] if exc["name"] == "Carbon dioxide, fossil"
    ):
        return dataset

    if "log parameters" not in dataset:
        dataset["log parameters"] = {}

    # subtract the biogenic CO2 amount to the initial fossil CO2 emission amount
    for exc in ws.biosphere(dataset, ws.equals("name", "Carbon dioxide, fossil")):
        dataset["log parameters"].update(
            {"initial amount of fossil CO2": exc["amount"]}
        )
        exc["amount"] -= amount_non_fossil_co2
        if exc["amount"] < 0:
            exc["amount"] = 0
        dataset["log parameters"].update({"new amount of fossil CO2": exc["amount"]})

    # add the non-fossil CO2 emission flow
    non_fossil_co2 = {
        "uncertainty type": 0,
        "amount": amount_non_fossil_co2,
        "type": "biosphere",
        "name": "Carbon dioxide, non-fossil",
        "unit": "kilogram",
        "categories": ("air",),
        "input": (
            "biosphere3",
            biosphere_flows[
                ("Carbon dioxide, non-fossil", "air", "unspecified", "kilogram")
            ],
        ),
    }

    dataset["log parameters"].update(
        {"new amount of biogenic CO2": amount_non_fossil_co2}
    )

    dataset["exchanges"].append(non_fossil_co2)

    return dataset


@lru_cache()
def add_boil_off_losses(vehicle, distance, loss_val):
    if vehicle == "truck":
        # average truck speed
        speed = 50
    else:
        # average ship speed
        speed = 36
    days = distance / speed / 24
    # boil-off losses, function of days in transit
    return np.power(1 + loss_val, days)


@lru_cache()
def add_pipeline_losses(distance, loss_val):
    # pipeline losses, function of distance
    return 1 + (loss_val * distance)


def add_other_losses(loss_val):
    return 1 + loss_val


def calculate_fuel_properties(amount, lhv, co2_factor, biogenic_share):
    """
    Calculate the fossil and non-fossil CO2 emissions and LHV for the given fuel
    properties and amount.
    """
    fossil_co2 = amount * lhv * co2_factor * (1 - biogenic_share)
    non_fossil_co2 = amount * lhv * co2_factor * biogenic_share
    weighted_lhv = amount * lhv
    return fossil_co2, non_fossil_co2, weighted_lhv


def update_dataset(dataset, supplier_key, amount):
    """
    Add a new exchange to the dataset for the given fuel and supplier, and update
    the LHV and CO2 fields.
    """
    exchange = {
        "uncertainty type": 0,
        "amount": amount,
        "product": supplier_key[2],
        "name": supplier_key[0],
        "unit": supplier_key[-1],
        "location": supplier_key[1],
        "type": "technosphere",
    }
    dataset["exchanges"].append(exchange)

    return dataset


def _update_fuels(scenario, version, system_model):
    fuels = Fuels(
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

    if any(
        x is not None
        for x in (
            scenario["iam data"].petrol_markets,
            scenario["iam data"].diesel_markets,
            scenario["iam data"].gas_markets,
            scenario["iam data"].hydrogen_markets,
        )
    ):
        fuels.generate_fuel_markets()
        fuels.adjust_fuel_conversion_efficiency()
        fuels.relink_datasets()
        scenario["database"] = fuels.database
        scenario["cache"] = fuels.cache
        scenario["index"] = fuels.index
    else:
        print("No fuel markets found in IAM data. Skipping.")

    return scenario


class Fuels(BaseTransformation):
    """
    Class that modifies fuel inventories and markets in ecoinvent based on IAM output data.
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
        # ecoinvent version
        self.version = version
        # dictionary of crops with specifications
        self.crops_props = get_crops_properties()
        # list to store markets that will be created
        self.new_fuel_markets = {}
        # dictionary to store mapping results, to avoid redundant effort
        self.cached_suppliers = {}
        self.biosphere_flows = get_biosphere_code(self.version)
        self.fuel_groups = fetch_mapping(FUEL_GROUPS)
        self.rev_fuel_groups = {
            sub_type: main_type
            for main_type, sub_types in self.fuel_groups.items()
            for sub_type in sub_types
        }

        self.iam_fuel_markets = self.iam_data.production_volumes.sel(
            variables=[
                g
                for g in [
                    item
                    for sublist in list(self.fuel_groups.values())
                    for item in sublist
                ]
                if g
                in self.iam_data.production_volumes.coords["variables"].values.tolist()
            ]
        )

        self.fuel_efficiencies = xr.DataArray(
            dims=["variables"], coords={"variables": []}
        )
        for efficiency in [
            self.iam_data.petrol_efficiencies,
            self.iam_data.diesel_efficiencies,
            self.iam_data.gas_efficiencies,
            self.iam_data.hydrogen_efficiencies,
        ]:
            if efficiency is not None:
                self.fuel_efficiencies = xr.concat(
                    [self.fuel_efficiencies, efficiency],
                    dim="variables",
                )

    def find_transport_activity(
        self, items_to_look_for: List[str], items_to_exclude: List[str], loc: str
    ) -> Tuple[str, str, str, str]:
        """Find the transport activity that is most similar to the given activity.
        This is done by looking for the most similar activity in the database.
        """

        try:
            dataset = ws.get_one(
                self.database,
                *[ws.contains("name", i) for i in items_to_look_for],
                ws.doesnt_contain_any("name", items_to_exclude),
                ws.equals("location", loc),
            )
        except ws.NoResults:
            dataset = ws.get_one(
                self.database,
                *[ws.contains("name", i) for i in items_to_look_for],
                ws.doesnt_contain_any("name", items_to_exclude),
            )

        return (
            dataset["name"],
            dataset["reference product"],
            dataset["unit"],
            dataset["location"],
        )

    def find_suppliers(
        self, name: str, ref_prod: str, unit: str, loc: str, exclude=None
    ) -> Dict[Tuple[Any, Any, Any, Any], float]:
        """
        Return a list of potential suppliers given a name, reference product,
        unit and location, with their respective supply share (based on production volumes).

        :param name: the name of the activity
        :param ref_prod: the reference product of the activity
        :param unit: the unit of the activity
        :param loc: the location of the activity
        :param exclude: a list of activities to exclude from the search
        :return: a dictionary of potential suppliers with their respective supply share
        """

        # if we find a result in the cache dictionary, return it
        if exclude is None:
            exclude = []
        key = (name, ref_prod, loc)
        if key in self.cached_suppliers:
            return self.cached_suppliers[key]

        ecoinvent_regions = self.geo.iam_to_ecoinvent_location(loc)
        # search first for suppliers in `loc`, but also in the ecoinvent
        # locations to are comprised in `loc`, and finally in `RER`, `RoW` and `GLO`,
        possible_locations = [
            [loc],
            ecoinvent_regions,
            ["RER"],
            ["RoW"],
            ["GLO"],
            ["CH"],
        ]
        suppliers, counter = [], 0

        # while we do not find a result
        while len(suppliers) == 0:
            try:
                suppliers = list(
                    get_suppliers_of_a_region(
                        database=self.database,
                        locations=possible_locations[counter],
                        names=[name] if isinstance(name, str) else name,
                        reference_prod=ref_prod,
                        unit=unit,
                        exclude=exclude,
                    )
                )
                counter += 1
            except IndexError as err:
                raise IndexError(
                    f"Could not find any supplier for {name} {ref_prod} in {possible_locations}."
                ) from err

        suppliers = [s for s in suppliers if s]  # filter out empty lists

        # find production volume-based share
        suppliers = get_shares_from_production_volume(suppliers)

        # store the result in cache for next time
        self.cached_suppliers[key] = suppliers

        return suppliers

    def generate_hydrogen_activities(self) -> None:
        """
        Defines regional variants for hydrogen production, but also different supply
        chain designs:
        * by truck (500 km), gaseous, liquid and LOHC
        * by reassigned CNG pipeline (500 km), gaseous, with and without inhibitors
        * by dedicated H2 pipeline (500 km), gaseous
        * by ship, liquid (2000 km)

        For truck and pipeline supply chains, we assume a transmission and a distribution part, for which
        we have specific pipeline designs. We also assume a means for regional storage in between (salt cavern).
        We apply distance-based losses along the way.

        Most of these supply chain design options are based on the work of:
        * Wulf C, Reuß M, Grube T, Zapp P, Robinius M, Hake JF, et al.
          Life Cycle Assessment of hydrogen transport and distribution options.
          J Clean Prod 2018;199:431–43. https://doi.org/10.1016/j.jclepro.2018.07.180.
        * Hank C, Sternberg A, Köppel N, Holst M, Smolinka T, Schaadt A, et al.
          Energy efficiency and economic assessment of imported energy carriers based on renewable electricity.
          Sustain Energy Fuels 2020;4:2256–73. https://doi.org/10.1039/d0se00067a.
        * Petitpas G. Boil-off losses along the LH2 pathway. US Dep Energy Off Sci Tech Inf 2018.


        """

        hydrogen_sources = fetch_mapping(HYDROGEN_SOURCES)

        for hydrogen_type, hydrogen_vars in hydrogen_sources.items():
            hydrogen_activity_name = hydrogen_sources[hydrogen_type].get("name")
            hydrogen_efficiency_variable = hydrogen_sources[hydrogen_type].get("var")
            hydrogen_feedstock_name = hydrogen_sources[hydrogen_type].get(
                "feedstock name"
            )
            hydrogen_feedstock_unit = hydrogen_sources[hydrogen_type].get(
                "feedstock unit"
            )
            efficiency_floor_value = hydrogen_sources[hydrogen_type].get("floor value")

            new_ds = self.fetch_proxies(
                name=hydrogen_activity_name,
                ref_prod="hydrogen",
                production_variable=hydrogen_efficiency_variable,
                exact_name_match=False,
            )

            for region, dataset in new_ds.items():
                # find current energy consumption in dataset
                initial_energy_consumption = sum(
                    exc["amount"]
                    for exc in dataset["exchanges"]
                    if exc["unit"] == hydrogen_feedstock_unit
                    and hydrogen_feedstock_name in exc["name"]
                    and exc["type"] != "production"
                )

                # add it to "log parameters"
                if "log parameters" not in dataset:
                    dataset["log parameters"] = {}

                dataset["log parameters"].update(
                    {
                        "initial energy input for hydrogen production": initial_energy_consumption
                    }
                )

                # Fetch the efficiency change of the
                # electrolysis process over time,
                # according to the IAM scenario,
                # if available.

                if (
                    hydrogen_efficiency_variable
                    in self.fuel_efficiencies.variables.values
                ):
                    # Find scaling factor compared to 2020
                    scaling_factor = 1 / self.find_iam_efficiency_change(
                        data=self.fuel_efficiencies,
                        variable=hydrogen_efficiency_variable,
                        location=region,
                    )

                    # new energy consumption
                    new_energy_consumption = scaling_factor * initial_energy_consumption

                    # set a floor value/kg H2
                    new_energy_consumption = max(
                        new_energy_consumption, efficiency_floor_value
                    )

                else:
                    if hydrogen_type == "from electrolysis":
                        # get the electricity consumption
                        new_energy_consumption = (
                            adjust_electrolysis_electricity_requirement(self.year)
                        )
                    else:
                        new_energy_consumption = None

                if new_energy_consumption is not None:
                    # recalculate scaling factor
                    scaling_factor = new_energy_consumption / initial_energy_consumption

                    if not np.isnan(scaling_factor) and scaling_factor > 0.0:
                        # rescale the fuel consumption exchange
                        dataset = rescale_exchanges(
                            dataset,
                            scaling_factor,
                            technosphere_filters=[
                                ws.contains("name", hydrogen_feedstock_name),
                                ws.equals("unit", hydrogen_feedstock_unit),
                            ],
                        )

                        # add it to "log parameters"
                        if "log parameters" not in dataset:
                            dataset["log parameters"] = {}

                        # add it to "log parameters"
                        dataset["log parameters"].update(
                            {
                                "new energy input for hydrogen production": new_energy_consumption
                            }
                        )

                    self.write_log(dataset)

                    string = f" The electricity input per kg of H2 has been adapted to the year {self.year}."
                    if "comment" in dataset:
                        dataset["comment"] += string
                    else:
                        dataset["comment"] = string

                    dataset["comment"] += (
                        "Region-specific hydrogen production dataset "
                        "generated by `premise`. "
                    )

            self.database.extend(new_ds.values())
            self.add_to_index(new_ds.values())

        # loss coefficients for hydrogen supply
        losses = fetch_mapping(HYDROGEN_SUPPLY_LOSSES)

        supply_chain_scenarios = fetch_mapping(SUPPLY_CHAIN_SCENARIOS)

        for act in [
            "hydrogen embrittlement inhibition",
            "geological hydrogen storage",
            "hydrogen refuelling station",
        ]:
            new_ds = self.fetch_proxies(name=act, ref_prod=" ")

            for k, dataset in new_ds.items():
                for exc in ws.production(dataset):
                    if "input" in exc:
                        del exc["input"]

                new_ds[k] = self.relink_technosphere_exchanges(
                    dataset,
                )

            self.database.extend(new_ds.values())

            # add to log
            for dataset in list(new_ds.values()):
                self.write_log(dataset)
                # add it to list of created datasets
                self.add_to_index(dataset)

        for region in self.regions:
            for hydrogen_type, hydrogen_vars in hydrogen_sources.items():
                for vehicle, config in supply_chain_scenarios.items():
                    for state in config["state"]:
                        for distance in config["distance"]:
                            # dataset creation
                            dataset: dict[
                                str,
                                Union[
                                    Union[
                                        str, list[dict[str, Union[int, str]]], ndarray
                                    ],
                                    Any,
                                ],
                            ] = {
                                "location": region,
                                "name": f"hydrogen supply, {hydrogen_type}, by {vehicle}, as {state}, over {distance} km",
                                "reference product": "hydrogen, 700 bar",
                                "unit": "kilogram",
                                "database": self.database[1]["database"],
                                "code": str(uuid.uuid4().hex),
                                "comment": "Dataset representing hydrogen supply, generated by `premise`.",
                                "exchanges": [
                                    {
                                        "uncertainty type": 0,
                                        "loc": 1,
                                        "amount": 1,
                                        "type": "production",
                                        "production volume": 1,
                                        "product": "hydrogen, 700 bar",
                                        "name": f"hydrogen supply, {hydrogen_type}, "
                                        f"by {vehicle}, as {state}, over {distance} km",
                                        "unit": "kilogram",
                                        "location": region,
                                    }
                                ],
                            }

                            # transport
                            dataset = self.add_hydrogen_transport(
                                dataset, config, region, distance, vehicle
                            )

                            # need for inhibitor and purification if CNG pipeline
                            # electricity for purification: 2.46 kWh/kg H2
                            if vehicle == "CNG pipeline":
                                dataset = self.add_hydrogen_inhibitor(dataset, region)

                            if "regional storage" in config:
                                dataset = self.add_hydrogen_regional_storage(
                                    dataset, region, config
                                )

                            # electricity for compression
                            if state in ["gaseous", "liquid"]:
                                dataset = self.add_compression_electricity(
                                    state, vehicle, distance, region, dataset
                                )

                            # electricity for hydrogenation, dehydrogenation and
                            # compression at delivery
                            if state == "liquid organic compound":
                                dataset = self.add_hydrogenation_energy(region, dataset)

                            dataset = self.add_hydrogen_input_and_losses(
                                hydrogen_vars,
                                region,
                                losses,
                                vehicle,
                                state,
                                distance,
                                dataset,
                            )

                            # add fuelling station, including storage tank
                            dataset["exchanges"].append(
                                self.add_h2_fuelling_station(region)
                            )

                            # add pre-cooling
                            dataset = self.add_pre_cooling_electricity(dataset, region)

                            dataset = self.relink_technosphere_exchanges(
                                dataset,
                            )

                            self.database.append(dataset)

                            # add to log
                            self.write_log(dataset)

                            # add it to list of created datasets
                            self.add_to_index(dataset)

    def add_hydrogen_transport(
        self,
        dataset: Dict[str, Any],
        config: Dict[str, Any],
        region: str,
        distance: float,
        vehicle: str,
    ) -> Dict[str, Any]:
        """
        Adds hydrogen transport exchanges to the given dataset.

        :param dataset: The dataset to modify.
        :param config: The configuration for the vehicle transport.
        :param region: The region of the dataset.
        :param distance: The distance traveled.
        :param vehicle: The type of vehicle used.
        :return: The modified dataset.
        """
        for transport in config["vehicle"]:
            transport_name = transport["name"]
            transport_ref_prod = transport["reference product"]
            transport_unit = transport["unit"]

            suppliers = self.find_suppliers(
                name=transport_name,
                ref_prod=transport_ref_prod,
                unit=transport_unit,
                loc=region,
            )

            for supplier, share in suppliers.items():
                if supplier[-1] == "ton kilometer":
                    amount = distance * share / 1000
                else:
                    amount = distance * share / 2 * (1 / eval(config["lifetime"]))

                exchange = {
                    "uncertainty type": 0,
                    "amount": amount,
                    "type": "technosphere",
                    "product": supplier[2],
                    "name": supplier[0],
                    "unit": supplier[-1],
                    "location": supplier[1],
                    "comment": f"Transport over {distance} km by {vehicle}. ",
                }

                dataset["exchanges"].append(exchange)

            comment = f"Transport over {distance} km by {vehicle}. "

            if "comment" in dataset:
                dataset["comment"] += comment
            else:
                dataset["comment"] = comment

        return dataset

    def add_hydrogen_input_and_losses(
        self, hydrogen_activity, region, losses, vehicle, state, distance, dataset
    ):
        # fetch the H2 production activity
        h2_ds = list(
            self.find_suppliers(
                name=hydrogen_activity["name"],
                ref_prod="hydrogen",
                unit="kilogram",
                loc=region,
            ).keys()
        )[0]

        # include losses along the way
        string = ""
        total_loss = 1
        for loss, val in losses[vehicle][state].items():
            val = float(val)

            if loss == "boil-off":
                total_loss *= add_boil_off_losses(vehicle, distance, val)
                string += f"Boil-off losses: {int((total_loss - 1) * 100)}%. "

            elif loss == "pipeline_leak":
                total_loss *= add_pipeline_losses(distance, val)
                string += f"Pipeline losses: {int((total_loss - 1) * 100)}%. "
            else:
                total_loss *= add_other_losses(val)
                string += f"{loss} losses: {int(val * 100)}%. "

        dataset["exchanges"].append(
            {
                "uncertainty type": 0,
                "amount": 1 * total_loss,
                "type": "technosphere",
                "product": h2_ds[2],
                "name": h2_ds[0],
                "unit": h2_ds[3],
                "location": region,
            }
        )

        # adds losses as hydrogen emission to air
        dataset["exchanges"].append(
            {
                "uncertainty type": 0,
                "amount": total_loss - 1,
                "type": "biosphere",
                "name": "Hydrogen",
                "unit": "kilogram",
                "categories": ("air",),
                "input": (
                    "biosphere3",
                    self.biosphere_flows[
                        (
                            "Hydrogen",
                            "air",
                            "unspecified",
                            "kilogram",
                        )
                    ],
                ),
            }
        )

        if "comment" in dataset:
            dataset["comment"] += string
        else:
            dataset["comment"] = string

        if "log parameters" not in dataset:
            dataset["log parameters"] = {}
        dataset["log parameters"].update(
            {"hydrogen distribution losses": total_loss - 1}
        )

        return dataset

    def add_hydrogenation_energy(
        self, region: str, dataset: Dict[str, Any]
    ) -> Dict[str, Any]:
        """
        Adds hydrogenation and dehydrogenation activities, as well as compression at delivery,
        to a dataset for a given region.

        :param region: The region for which to add the activities.
        :param dataset: The dataset to modify.
        :return: The modified dataset.

        :raises ValueError: If no hydrogenation activity is found for the specified region.

        """

        try:
            # fetch the H2 hydrogenation activity
            hydrogenation_ds = list(
                self.find_suppliers(
                    name="hydrogenation of hydrogen",
                    ref_prod="hydrogenation",
                    unit="kilogram",
                    loc=region,
                ).keys()
            )[0]

            # fetch the H2 de-hydrogenation activity
            dehydrogenation_ds = list(
                self.find_suppliers(
                    name="dehydrogenation of hydrogen",
                    ref_prod="dehydrogenation",
                    unit="kilogram",
                    loc=region,
                ).keys()
            )[0]

        except ws.NoResults as err_noresults:
            raise ValueError(
                f"No hydrogenation activity found for region {region}"
            ) from err_noresults
        except ws.MultipleResults as err_multipleresults:
            raise ValueError(
                f"Multiple hydrogenation activities found for region {region}"
            ) from err_multipleresults

        dataset["exchanges"].extend(
            [
                {
                    "uncertainty type": 0,
                    "amount": 1,
                    "type": "technosphere",
                    "product": hydrogenation_ds[2],
                    "name": hydrogenation_ds[0],
                    "unit": hydrogenation_ds[3],
                    "location": region,
                },
                {
                    "uncertainty type": 0,
                    "amount": 1,
                    "type": "technosphere",
                    "product": dehydrogenation_ds[2],
                    "name": dehydrogenation_ds[0],
                    "unit": dehydrogenation_ds[3],
                    "location": region,
                },
            ]
        )

        # After dehydrogenation at ambient temperature at delivery
        # the hydrogen needs to be compressed up to 900 bar to be dispensed
        # in 700 bar storage tanks

        electricity_comp = get_compression_effort(25, 900, 1000)

        electricity_suppliers = self.find_suppliers(
            name="market group for electricity, low voltage",
            ref_prod="electricity, low voltage",
            unit="kilowatt hour",
            loc=region,
            exclude=["period"],
        )

        dataset["exchanges"].extend(
            {
                "uncertainty type": 0,
                "amount": electricity_comp * share,
                "type": "technosphere",
                "product": supplier[2],
                "name": supplier[0],
                "unit": supplier[-1],
                "location": supplier[1],
            }
            for supplier, share in electricity_suppliers.items()
        )

        string = (
            " Hydrogenation and dehydrogenation of hydrogen included. "
            "Compression at delivery after dehydrogenation also included."
        )
        if "comment" in dataset:
            dataset["comment"] += string
        else:
            dataset["comment"] = string

        if "log parameters" not in dataset:
            dataset["log parameters"] = {}
        dataset["log parameters"].update(
            {
                "electricity for hydrogen compression after dehydrogenation": electricity_comp
            }
        )

        return dataset

    def add_hydrogen_regional_storage(
        self, dataset: dict, region: str, config: dict
    ) -> dict:
        """

        Add a geological storage activity to the dataset for a given region.

        :param dataset: The dataset to modify.
        :param region: The region for which to add the activity.
        :param config: The configuration file for the analysis.
        :return: The modified dataset.

        """

        storage_ds = list(
            self.find_suppliers(
                name=config["regional storage"]["name"],
                ref_prod=config["regional storage"]["reference product"],
                unit=config["regional storage"]["unit"],
                loc=region,
            ).keys()
        )[0]

        dataset["exchanges"].append(
            {
                "uncertainty type": 0,
                "amount": 1,
                "type": "technosphere",
                "product": storage_ds[2],
                "name": storage_ds[0],
                "unit": storage_ds[3],
                "location": region,
                "comment": "Geological storage (salt cavern).",
            }
        )

        string = (
            " Geological storage is added. It includes 0.344 kWh for "
            "the injection and pumping of 1 kg of H2."
        )
        if "comment" in dataset:
            dataset["comment"] += string
        else:
            dataset["comment"] = string

        return dataset

    def add_hydrogen_inhibitor(self, dataset: dict, region: str) -> dict:
        """
        Adds hydrogen embrittlement inhibitor to the dataset for a given region.

        :param dataset: The dataset to modify.
        :param region: The region for which to add the activity.
        :return: The modified dataset.
        """

        inhibbitor_ds = list(
            self.find_suppliers(
                name="hydrogen embrittlement inhibition",
                ref_prod="hydrogen",
                unit="kilogram",
                loc=region,
            ).keys()
        )[0]

        dataset["exchanges"].append(
            {
                "uncertainty type": 0,
                "amount": 1,
                "type": "technosphere",
                "product": inhibbitor_ds[2],
                "name": inhibbitor_ds[0],
                "unit": inhibbitor_ds[3],
                "location": region,
                "comment": "Injection of an inhibiting gas (oxygen) "
                "to prevent embrittlement of metal. ",
            }
        )

        string = (
            "2.46 kWh/kg H2 is needed to purify the hydrogen from the inhibiting gas. "
            "The recovery rate for hydrogen after separation from the inhibitor gas is 93%. "
        )
        if "comment" in dataset:
            dataset["comment"] += string
        else:
            dataset["comment"] = string

        return dataset

    def add_compression_electricity(
        self, state: str, vehicle: str, distance: float, region: str, dataset: dict
    ) -> dict:
        """
        Add the electricity needed for the compression of hydrogen.

        :param state: The state of the hydrogen (gaseous or liquid).
        :param vehicle: The vehicle used for transport (truck or pipeline).
        :param distance: The distance travelled by the vehicle.
        :param region: The region for which to add the activity.
        :param dataset: The dataset to modify.
        :return: The modified dataset.

        """

        # if gaseous
        # if transport by truck, compression from 25 bar to 500 bar for the transport
        # and from 500 bar to 900 bar for dispensing in 700 bar storage tanks

        # if transport by pipeline, initial compression from 25 bar to 100 bar
        # and 0.6 kWh re-compression every 250 km
        # and finally from 100 bar to 900 bar for dispensing in 700 bar storage tanks

        # if liquid
        # liquefaction electricity need
        # currently, 12 kWh/kg H2
        # midterm, 8 kWh/ kg H2
        # by 2050, 6 kWh/kg H2

        if state == "gaseous":
            if vehicle == "truck":
                electricity_comp = get_compression_effort(25, 500, 1000)
                electricity_comp += get_compression_effort(500, 900, 1000)
            else:
                electricity_comp = get_compression_effort(25, 100, 1000) + (
                    0.6 * distance / 250
                )
                electricity_comp += get_compression_effort(100, 900, 1000)

            string = (
                f" {electricity_comp} kWh is added to compress from 25 bar 100 bar (if pipeline)"
                f"or 500 bar (if truck), and then to 900 bar to dispense in storage tanks at 700 bar. "
                " Additionally, if transported by pipeline, there is re-compression (0.6 kWh) every 250 km."
            )

        else:
            electricity_comp = np.clip(
                np.interp(
                    self.year,
                    [2020, 2035, 2050],
                    [12, 8, 6],
                ),
                12,
                6,
            )

            string = f" {electricity_comp} kWh is added to liquefy the hydrogen. "

        suppliers = self.find_suppliers(
            name="market group for electricity, low voltage",
            ref_prod="electricity, low voltage",
            unit="kilowatt hour",
            loc=region,
            exclude=["period"],
        )

        new_exc = []
        for supplier, share in suppliers.items():
            new_exc.append(
                {
                    "uncertainty type": 0,
                    "amount": electricity_comp * share,
                    "type": "technosphere",
                    "product": supplier[2],
                    "name": supplier[0],
                    "unit": supplier[-1],
                    "location": supplier[1],
                }
            )

        dataset["exchanges"].extend(new_exc)

        if "comment" in dataset:
            dataset["comment"] += string
        else:
            dataset["comment"] = string

        if "log parameters" not in dataset:
            dataset["log parameters"] = {}
        dataset["log parameters"].update(
            {"electricity for hydrogen compression": electricity_comp}
        )

        return dataset

    @lru_cache()
    def add_h2_fuelling_station(self, region: str) -> dict:
        """
        Add the hydrogen fuelling station.

        :param region: The region for which to add the activity.
        :return: The modified dataset.

        """

        ds_h2_station = list(
            self.find_suppliers(
                name="hydrogen refuelling station",
                ref_prod="hydrogen",
                unit="unit",
                loc=region,
            ).keys()
        )[0]

        return {
            "uncertainty type": 0,
            "amount": 1 / (600 * 365 * 40),  # 1 over lifetime: 40 years, 600 kg H2/day
            "type": "technosphere",
            "product": ds_h2_station[2],
            "name": ds_h2_station[0],
            "unit": ds_h2_station[3],
            "location": region,
        }

    def add_pre_cooling_electricity(self, dataset: dict, region: str) -> dict:
        """
        Add the electricity needed for pre-cooling the hydrogen.

        :param dataset: The dataset to modify.
        :param region: The region for which to add the activity.
        :return: The modified dataset.
        """

        # finally, add pre-cooling
        # is needed before filling vehicle tanks
        # as the hydrogen is pumped, the ambient temperature
        # vaporizes the gas, and because of the Thomson-Joule effect,
        # the gas temperature increases.
        # Hence, a refrigerant is needed to keep the H2 as low as
        # -30 C during pumping.

        # https://www.osti.gov/servlets/purl/1422579 gives us a formula
        # to estimate pre-cooling electricity need
        # it requires a capacity utilization for the fuelling station
        # as well as an ambient temperature
        # we will use a temp of 25 C
        # and a capacity utilization going from 10 kg H2/day in 2020
        # to 150 kg H2/day in 2050

        t_amb = 25
        cap_util = np.interp(self.year, [2020, 2050, 2100], [10, 150, 150])
        el_pre_cooling = get_pre_cooling_energy(t_amb, float(cap_util))

        suppliers = self.find_suppliers(
            name="market group for electricity, low voltage",
            ref_prod="electricity, low voltage",
            unit="kilowatt hour",
            loc=region,
            exclude=["period"],
        )

        for supplier, share in suppliers.items():
            dataset["exchanges"].append(
                {
                    "uncertainty type": 0,
                    "amount": el_pre_cooling * share,
                    "type": "technosphere",
                    "product": supplier[2],
                    "name": supplier[0],
                    "unit": supplier[-1],
                    "location": supplier[1],
                }
            )

        string = (
            f"Pre-cooling electricity is considered ({el_pre_cooling}), "
            f"assuming an ambiant temperature of {t_amb}C "
            f"and a capacity utilization for the fuel station of {cap_util} kg/day."
        )
        if "comment" in dataset:
            dataset["comment"] += string
        else:
            dataset["comment"] = string

        if "log parameters" not in dataset:
            dataset["log parameters"] = {}
        dataset["log parameters"].update(
            {"electricity for hydrogen pre-cooling": el_pre_cooling}
        )

        return dataset

    def generate_biogas_activities(self):
        """
        Generate biogas activities.
        """

        fuel_activities = fetch_mapping(METHANE_SOURCES)

        for fuel, activities in fuel_activities.items():
            for activity in activities:
                if fuel == "methane, synthetic":
                    original_ds = self.fetch_proxies(
                        name=activity,
                        ref_prod=" ",
                        delete_original_dataset=False,
                        empty_original_activity=False,
                        exact_name_match=False,
                        relink=False,
                    )

                    for co2_type in [
                        (
                            "carbon dioxide, captured from atmosphere, with a sorbent-based direct air capture system, 100ktCO2, with waste heat, and grid electricity",
                            "carbon dioxide, captured from atmosphere",
                            "waste heat",
                        ),
                        (
                            "carbon dioxide, captured from atmosphere, with a solvent-based direct air capture system, 1MtCO2, with heat pump heat, and grid electricity",
                            "carbon dioxide, captured from atmosphere",
                            "heat pump heat",
                        ),
                    ]:
                        new_ds = copy.deepcopy(original_ds)

                        for region, dataset in new_ds.items():
                            dataset["code"] = str(uuid.uuid4().hex)
                            dataset["name"] += f", using {co2_type[2]}"
                            for prod in ws.production(dataset):
                                prod["name"] = dataset["name"]

                                if "input" in prod:
                                    del prod["input"]

                            for exc in ws.technosphere(dataset):
                                if (
                                    "carbon dioxide, captured from atmosphere"
                                    in exc["name"].lower()
                                ):
                                    # store amount
                                    co2_amount = exc["amount"]

                                    try:
                                        # add new exchanges
                                        dac_suppliers = self.find_suppliers(
                                            name=co2_type[0],
                                            ref_prod=co2_type[1],
                                            unit="kilogram",
                                            loc=region,
                                        )

                                    except IndexError:
                                        dac_suppliers = None

                                    if dac_suppliers:
                                        # remove exchange
                                        dataset["exchanges"].remove(exc)

                                        dataset["exchanges"].extend(
                                            {
                                                "uncertainty type": 0,
                                                "amount": co2_amount * share,
                                                "type": "technosphere",
                                                "product": supplier[2],
                                                "name": supplier[0],
                                                "unit": supplier[-1],
                                                "location": supplier[1],
                                            }
                                            for supplier, share in dac_suppliers.items()
                                        )

                            for exc in ws.technosphere(dataset):
                                if (
                                    "methane, from electrochemical methanation"
                                    in exc["name"]
                                ):
                                    exc["name"] += f", using {co2_type[2]}"
                                    exc["location"] = dataset["location"]

                                    dataset["name"] = dataset["name"].replace(
                                        "from electrochemical methanation",
                                        f"from electrochemical methanation "
                                        f"(H2 from electrolysis, CO2 from DAC "
                                        f"using {co2_type[2]})",
                                    )

                                    for prod in ws.production(dataset):
                                        prod["name"] = prod["name"].replace(
                                            "from electrochemical methanation",
                                            f"from electrochemical methanation "
                                            f"(H2 from electrolysis, CO2 from DAC "
                                            f"using {co2_type[2]})",
                                        )

                        self.database.extend(new_ds.values())

                        # add to log
                        for new_dataset in list(new_ds.values()):
                            self.write_log(new_dataset)

                            # add it to list of created datasets
                            self.add_to_index(new_dataset)
                else:
                    original_ds = self.fetch_proxies(
                        name=activity, ref_prod=" ", exact_name_match=False
                    )
                    new_ds = copy.deepcopy(original_ds)

                    for region, dataset in new_ds.items():
                        dataset["code"] = str(uuid.uuid4().hex)
                        for exc in ws.production(dataset):
                            if "input" in exc:
                                exc.pop("input")

                        new_ds[region] = self.relink_technosphere_exchanges(
                            dataset,
                        )

                    self.database.extend(new_ds.values())

                    # add to log
                    for new_dataset in list(new_ds.values()):
                        self.write_log(new_dataset)
                        # add it to list of created datasets
                        self.add_to_index(new_dataset)

    def generate_synthetic_fuel_activities(self):
        """
        Generate synthetic fuel activities.
        """

        fuel_activities = fetch_mapping(LIQUID_FUEL_SOURCES)

        for activities in fuel_activities.values():
            for activity in activities:
                new_ds = self.fetch_proxies(
                    name=activity,
                    ref_prod=" ",
                    relink=False,
                )
                for region, dataset in new_ds.items():
                    for exc in ws.production(dataset):
                        if "input" in exc:
                            del exc["input"]

                    for exc in ws.technosphere(dataset):
                        if "carbon dioxide, captured from atmosphere" in exc["name"]:
                            # store amount
                            co2_amount = exc["amount"]

                            try:
                                # add new exchanges
                                dac_suppliers = self.find_suppliers(
                                    name="carbon dioxide, captured from atmosphere, with a solvent-based direct air capture "
                                    "system, 1MtCO2, with heat pump heat, and grid electricity",
                                    ref_prod="carbon dioxide, captured from atmosphere",
                                    unit="kilogram",
                                    loc=region,
                                )

                            except IndexError:
                                dac_suppliers = None

                            if dac_suppliers:
                                # remove exchange
                                dataset["exchanges"].remove(exc)

                                dataset["exchanges"].extend(
                                    {
                                        "uncertainty type": 0,
                                        "amount": co2_amount * share,
                                        "type": "technosphere",
                                        "product": supplier[2],
                                        "name": supplier[0],
                                        "unit": supplier[-1],
                                        "location": supplier[1],
                                    }
                                    for supplier, share in dac_suppliers.items()
                                )

                    dataset = self.relink_technosphere_exchanges(
                        dataset,
                    )

                    self.database.append(dataset)

                    # add to log
                    self.write_log(dataset)
                    # add it to list of created datasets
                    self.add_to_index(dataset)

    def adjust_land_use(self, dataset: dict, region: str, crop_type: str) -> dict:
        """
        Adjust land use.

        :param dataset: dataset to adjust
        :param region: region of the dataset
        :param crop_type: crop type of the dataset
        :return: adjusted dataset

        """
        string = ""
        land_use = 0

        for exc in dataset["exchanges"]:
            # we adjust the land use
            if exc["type"] == "biosphere" and exc["name"].startswith("Occupation"):
                if "LHV [MJ/kg as received]" in dataset:
                    lower_heating_value = dataset["LHV [MJ/kg as received]"]
                else:
                    lower_heating_value = dataset.get("LHV [MJ/kg dry]", 0)

                # Ha/GJ
                land_use = (
                    self.iam_data.land_use.sel(region=region, variables=crop_type)
                    .interp(year=self.year)
                    .values
                )

                # replace NA values with 0
                if np.isnan(land_use):
                    land_use = 0

                if land_use > 0:
                    # HA to m2
                    land_use *= 10000
                    # m2/GJ to m2/MJ
                    land_use /= 1000
                    # m2/kg, as received
                    land_use *= lower_heating_value
                    # update exchange value
                    exc["amount"] = float(land_use)

                    string = (
                        f"The land area occupied has been modified to {land_use}, "
                        f"to be in line with the scenario {self.scenario} of {self.model.upper()} "
                        f"in {self.year} in the region {region}. "
                    )

        if string and land_use:
            if "comment" in dataset:
                dataset["comment"] += string
            else:
                dataset["comment"] = string

            if "log parameters" not in dataset:
                dataset["log parameters"] = {}
            dataset["log parameters"].update(
                {
                    "land footprint": land_use,
                }
            )

        return dataset

    def adjust_land_use_change_emissions(
        self,
        dataset: dict,
        region: str,
        crop_type: str,
    ) -> dict:
        """
        Adjust land use change emissions to crop farming dataset
        if the variable is provided by the IAM.

        :param dataset: dataset to adjust
        :param region: region of the dataset
        :param crop_type: crop type of the dataset
        :return: adjusted dataset

        """

        # then, we should include the Land Use Change-induced CO2 emissions
        # those are given in kg CO2-eq./GJ of primary crop energy

        # kg CO2/GJ
        land_use_co2 = (
            self.iam_data.land_use_change.sel(region=region, variables=crop_type)
            .interp(year=self.year)
            .values
        )

        # replace NA values with 0
        if np.isnan(land_use_co2):
            land_use_co2 = 0

        if land_use_co2 > 0:
            # lower heating value, as received
            if "LHV [MJ/kg as received]" in dataset:
                lower_heating_value = dataset["LHV [MJ/kg as received]"]
            else:
                lower_heating_value = dataset.get("LHV [MJ/kg dry]", 0)

            # kg CO2/MJ
            land_use_co2 /= 1000
            land_use_co2 *= lower_heating_value

            land_use_co2_exc = {
                "uncertainty type": 0,
                "loc": float(land_use_co2),
                "amount": float(land_use_co2),
                "type": "biosphere",
                "name": "Carbon dioxide, from soil or biomass stock",
                "unit": "kilogram",
                "input": (
                    "biosphere3",
                    self.biosphere_flows[
                        (
                            "Carbon dioxide, from soil or biomass stock",
                            "air",
                            "non-urban air or from high stacks",
                            "kilogram",
                        )
                    ],
                ),
                "categories": (
                    "air",
                    "non-urban air or from high stacks",
                ),
            }
            dataset["exchanges"].append(land_use_co2_exc)

            string = (
                f"{land_use_co2} kg of land use-induced CO2 has been added by premise, "
                f"to be in line with the scenario {self.scenario} of {self.model.upper()} "
                f"in {self.year} in the region {region}."
            )

            if "comment" in dataset:
                dataset["comment"] += string
            else:
                dataset["comment"] = string

            if "log parameters" not in dataset:
                dataset["log parameters"] = {}
            dataset["log parameters"].update(
                {
                    "land use CO2": land_use_co2,
                }
            )

        return dataset

    def adjust_fuel_conversion_efficiency(self):
        """
        Adjust the input to output fuel conversion efficiency.
        """

        for fuel, activities in self.fuel_map.items():
            for activity in activities:
                for ds in ws.get_many(self.database, ws.equals("name", activity)):
                    variable = self.rev_fuel_map.get(activity)
                    scaling_factor = 1.0
                    if variable in self.fuel_efficiencies.coords["variables"]:
                        if ds["location"] in self.regions:
                            region = ds["location"]
                        else:
                            region = self.ecoinvent_to_iam_loc[ds["location"]]

                        scaling_factor = (
                            self.fuel_efficiencies.sel(
                                variables=variable,
                                region=region,
                            )
                            .interp(year=self.year)
                            .values
                        )
                    if (
                        scaling_factor != 1.0
                        and "market for" not in ds["name"]
                        and "fuel conversion efficiency"
                        not in ds.get("log parameters", {})
                    ):
                        rescale_exchanges(
                            ds,
                            1 / scaling_factor,
                        )

                        if "log parameters" not in ds:
                            ds["log parameters"] = {}

                        ds["log parameters"].update(
                            {
                                "fuel conversion efficiency": 1 / scaling_factor,
                            }
                        )

                        # update comment section
                        txt = (
                            f" The inputs of this dataset have been multiplied by {1 / scaling_factor}"
                            f"to reflect changes in efficiency according to {self.model.upper()} - {self.scenario.upper()}."
                        )
                        if "comment" in ds:
                            ds["comment"] += txt
                        else:
                            ds["comment"] = txt
                        self.write_log(ds, status="updated")

    def get_production_label(self, crop_type: str) -> [str, None]:
        """
        Get the production label for the dataset.
        """
        try:
            return [
                i
                for i in self.iam_fuel_markets.coords["variables"].values.tolist()
                if crop_type.lower() in i.lower()
            ][0]
        except IndexError:
            return None

    def should_adjust_land_use(self, dataset: dict, crop_type: str) -> bool:
        """
        Check if the dataset should be adjusted for land use.
        """

        if self.iam_data.land_use is None:
            return False
        return (
            any(i in dataset["name"].lower() for i in ("farming and supply",))
            and crop_type.lower() in self.iam_data.land_use.variables.values
            and not any(
                i in dataset["name"].lower() for i in ["straw", "residue", "stover"]
            )
        )

    def should_adjust_land_use_change_emissions(
        self, dataset: dict, crop_type: str
    ) -> bool:
        """
        Check if the dataset should be adjusted for land use change emissions.
        """
        if self.iam_data.land_use_change is None:
            return False
        return (
            any(i in dataset["name"].lower() for i in ("farming and supply",))
            and crop_type.lower() in self.iam_data.land_use_change.variables.values
            and not any(
                i in dataset["name"].lower() for i in ["straw", "residue", "stover"]
            )
        )

    def generate_biofuel_activities(self):
        """
        Create region-specific biofuel datasets.
        Update the conversion efficiency.
        :return:
        """

        # Map regions to their respective climate types
        region_to_climate = fetch_mapping(REGION_CLIMATE_MAP)[self.model]

        # Map climate types to their respective crop types and crops
        crop_types = list(self.crops_props.keys())
        climates = set(region_to_climate.values())
        climate_to_crop_type = {
            clim: {
                crop_type: self.crops_props[crop_type]["crop_type"][self.model][clim]
                for crop_type in crop_types
            }
            for clim in climates
        }

        biofuel_activities = fetch_mapping(BIOFUEL_SOURCES)

        # List to store processed crops
        processed_crops = []

        for climate in ["tropical", "temperate"]:
            regions = [k for k, v in region_to_climate.items() if v == climate]
            for crop_type in climate_to_crop_type[climate]:
                specific_crop = climate_to_crop_type[climate][crop_type]
                if specific_crop in processed_crops:
                    continue
                processed_crops.append(specific_crop)

                # Skip processing corn if it's already been processed
                if specific_crop == "corn":
                    regions = list(region_to_climate.keys())

                # Get list of activities for the crop and biofuel type
                activities = biofuel_activities[crop_type][specific_crop]

                for activity in activities:
                    # Fetch dataset for activity and regions
                    new_datasets = self.fetch_proxies(
                        name=activity,
                        ref_prod=" ",
                        production_variable=self.get_production_label(
                            crop_type=crop_type
                        ),
                        regions=regions,
                    )

                    # Adjust land use for farming activities
                    new_datasets = {
                        region: (
                            self.adjust_land_use(ds, region, crop_type)
                            if self.should_adjust_land_use(ds, crop_type)
                            else ds
                        )
                        for region, ds in new_datasets.items()
                    }

                    # Adjust land use change emissions for farming activities
                    new_datasets = {
                        region: (
                            self.adjust_land_use_change_emissions(ds, region, crop_type)
                            if self.should_adjust_land_use_change_emissions(
                                ds, crop_type
                            )
                            else ds
                        )
                        for region, ds in new_datasets.items()
                    }

                    self.database.extend(new_datasets.values())

                    # add to log
                    for dataset in list(new_datasets.values()):
                        self.write_log(dataset)

                        # add it to list of created datasets
                        self.add_to_index(dataset)

    def get_fuel_mapping(self) -> dict:
        """
        Define filter functions that decide which wurst datasets to modify.
        :return: dictionary that contains filters and functions
        :rtype: dict
        """

        return {
            fuel: {
                "find_share": self.fetch_fuel_share,
                "fuel filters": self.fuel_map[fuel],
            }
            for fuel in self.iam_fuel_markets.variables.values
        }

    @lru_cache()
    def fetch_fuel_share(
        self, fuel: str, relevant_fuel_types: Tuple[str], region: str, period: int
    ) -> float:
        """
        Return the percentage of a specific fuel type in the fuel mix for a specific region.
        :param fuel: the name of the fuel to fetch the percentage for
        :param relevant_fuel_types: a list of relevant fuel types to include in the calculation
        :param region: the IAM region to fetch the data for
        :param period: the period to fetch the data for
        :return: the percentage of the specified fuel type in the fuel mix for the region
        """

        relevant_variables = [
            v
            for v in self.iam_fuel_markets.variables.values
            if any(x.lower() in v.lower() for x in relevant_fuel_types)
        ]

        fuel_share = (
            (
                self.iam_fuel_markets.sel(region=region, variables=fuel)
                / self.iam_fuel_markets.sel(
                    region=region, variables=relevant_variables
                ).sum(dim="variables")
            )
            .fillna(0)
            .interp(
                year=np.arange(self.year, self.year + period + 1),
                kwargs={"fill_value": "extrapolate"},
            )
            .mean(dim="year")
            .values
        )

        if np.isnan(fuel_share):
            print(
                f"Warning: incorrect fuel share for {fuel} in {region} (-> set to 0%)."
            )
            fuel_share = 0

        return float(fuel_share)

    def relink_activities_to_new_markets(self):
        """
        Links fuel input exchanges to new datasets
        with the appropriate IAM location.

        Does not return anything.
        """

        # Create set of activities that consume fuels
        created_markets = list(set(x[0] for x in self.new_fuel_markets))
        # add "market group for" activities
        created_markets.extend(
            [
                "market group for " + x.replace("market for ", "")
                for x in list(set(x[0] for x in self.new_fuel_markets))
            ]
        )
        # add import datasets
        created_markets.append("import from")
        # add natural gas, petrol and diesel production datasets
        created_markets.extend(
            [
                "natural gas production",
                "petrol production",
                "diesel production",
            ]
        )
        # transport datasets
        created_markets.append("transport, pipeline")
        # and pyrolysis-based H2
        created_markets.append("from methane pyrolysis")

        new_keys = {}
        for key, value in self.new_fuel_markets.items():
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

        self.new_fuel_markets.update(new_keys)

        # Get fuel markets and amounts
        fuel_markets = fetch_mapping(FUEL_MARKETS)

        list_items_to_ignore = [
            "blending",
            "lubricating oil production",
            "petrol production",
            "natural gas production",
        ]

        old_fuel_inputs = ["market for " + x for x in list(fuel_markets.keys())] + [
            "market for petrol, unleaded",
            "market for diesel",
        ]
        old_fuel_inputs.extend(
            ["market group for " + x for x in list(fuel_markets.keys())]
        )
        old_fuel_inputs.extend(
            [
                "market group for petrol, unleaded",
                "market group for diesel",
                "market group for diesel, low-sulfur",
            ]
        )
        old_fuel_inputs = list(set(old_fuel_inputs))

        # Iterate over datasets and update exchanges as necessary
        for dataset in ws.get_many(
            self.database,
            ws.exclude(ws.either(*[ws.contains("name", x) for x in created_markets])),
        ):
            # Check that a fuel input exchange is present
            # in the list of inputs
            # Check also for "market group for" inputs
            exchanges = list(
                ws.technosphere(
                    dataset,
                    ws.either(*[ws.contains("name", x) for x in old_fuel_inputs]),
                )
            )

            if len(exchanges) > 0:
                supplier_loc = (
                    dataset["location"]
                    if dataset["location"] in self.regions
                    else self.geo.ecoinvent_to_iam_location(dataset["location"])
                )

                amount_non_fossil_co2 = sum(
                    a["amount"]
                    * self.new_fuel_markets.get(
                        (
                            a["name"],
                            supplier_loc,
                        ),
                        {},
                    ).get("non-fossil CO2", 0)
                    for a in exchanges
                )

                if amount_non_fossil_co2 > 0 and not any(
                    x in dataset["name"].lower() for x in list_items_to_ignore
                ):
                    update_co2_emissions(
                        dataset, amount_non_fossil_co2, self.biosphere_flows
                    )
                    self.write_log(dataset, status="updated")

    def generate_fuel_supply_chains(self):
        """Duplicate fuel chains and make them IAM region-specific"""

        # hydrogen
        self.generate_hydrogen_activities()

        # biogas
        self.generate_biogas_activities()

        # synthetic fuels
        self.generate_synthetic_fuel_activities()

        # biofuels
        self.generate_biofuel_activities()

    def generate_world_fuel_market(
        self, dataset: dict, d_act: dict, prod_vars: list, period: int
    ) -> dict:
        """
        Generate the world fuel market for a given dataset and product variables.

        :param dataset: The dataset for which to generate the world fuel market.
        :param d_act: A dictionary of activity datasets, keyed by region.
        :param prod_vars: A list of product variables.
        :return: A tuple containing the final LHV, fossil CO2, and biogenic CO2 emissions for the world fuel market,


        This function generates the world fuel market exchanges for a given dataset and set of product variables.
        It first filters out non-production exchanges from the dataset, and then calculates the total production
        volume for the world using the given product variables. For each region, it calculates the share of the
        production volume and adds a technosphere exchange to the dataset with the appropriate share. It also
        calculates the total LHV, fossil CO2, and biogenic CO2 emissions for each region. Finally, it returns a
        tuple with the final LHV, fossil CO2, and biogenic CO2 emissions for the world fuel market, as well as the
        updated dataset with the world fuel market exchanges.

        """

        if period != 0:
            # this dataset is for a period of time
            dataset["name"] += f", {period}-year period"
            dataset["comment"] += (
                f" Average fuel mix over a {period}"
                f"-year period {self.year}-{self.year + period}."
            )
            for exc in ws.production(dataset):
                exc["name"] += f", {period}-year period"

        # Filter out non-production exchanges
        dataset["exchanges"] = [
            e for e in dataset["exchanges"] if e["type"] == "production"
        ]

        final_lhv, final_fossil_co2, final_biogenic_co2 = 0, 0, 0

        # Calculate share of production volume for each region
        for r in d_act.keys():
            if r == "World" or (dataset["name"], r) not in self.new_fuel_markets:
                continue

            share = (
                (
                    self.iam_fuel_markets.sel(region=r, variables=prod_vars).sum(
                        dim="variables"
                    )
                    / self.iam_fuel_markets.sel(
                        variables=prod_vars,
                        region=[
                            x
                            for x in self.iam_fuel_markets.region.values
                            if x != "World"
                        ],
                    ).sum(dim=["variables", "region"])
                )
                .interp(
                    year=np.arange(self.year, self.year + period + 1),
                    kwargs={"fill_value": "extrapolate"},
                )
                .mean(dim="year")
                .values
            )

            if np.isnan(share):
                print("Incorrect market share for", dataset["name"], "in", r)

            # Calculate total LHV, fossil CO2, and biogenic CO2 for the region
            fuel_market_key = (dataset["name"], r)

            # if key absent from self.new_fuel_markets, then it does not exist

            if fuel_market_key in self.new_fuel_markets and share > 0:
                # Add exchange for the region
                exchange = {
                    "uncertainty type": 0,
                    "amount": share,
                    "type": "technosphere",
                    "product": dataset["reference product"],
                    "name": dataset["name"],
                    "unit": dataset["unit"],
                    "location": r,
                }
                dataset["exchanges"].append(exchange)

                lhv = self.new_fuel_markets[fuel_market_key]["LHV"]
                co2_factor = self.new_fuel_markets[fuel_market_key]["fossil CO2"]
                biogenic_co2_factor = self.new_fuel_markets[fuel_market_key][
                    "non-fossil CO2"
                ]
                final_lhv += share * lhv
                final_fossil_co2 += share * co2_factor
                final_biogenic_co2 += share * biogenic_co2_factor

                dataset["log parameters"] = {}
                dataset["log parameters"]["fossil CO2 per kg fuel"] = final_fossil_co2
                dataset["log parameters"][
                    "non-fossil CO2 per kg fuel"
                ] = final_biogenic_co2
                dataset["log parameters"]["lower heating value"] = final_lhv

        return dataset

    def generate_regional_fuel_market(
        self,
        dataset: dict,
        fuel_providers: dict,
        prod_vars: list,
        vars_map: dict,
        fuel_category: str,
        region: str,
        activity: dict,
        period: int,
    ) -> dict:
        """
        Generate regional fuel market for a given dataset and fuel providers.

        :param dataset: The dataset for which to generate the regional fuel market.
        :param fuel_providers: A dictionary of fuel providers, keyed by product variable.
        :param prod_vars: A list of product variables.
        :param vars_map: A dictionary mapping product variables to fuel names.
        :param fuel_category: The fuel name.
        :param region: The region for which to generate the regional fuel market.
        :param activity: The activity dataset for the region.
        :return: A tuple containing the final LHV, fossil CO2, and biogenic CO2 emissions for the regional fuel market,
        as well as the updated dataset with the regional fuel market exchanges.

        """

        # Initialize variables
        fossil_co2, non_fossil_co2, final_lhv = [0, 0, 0]

        if period != 0:
            # this dataset is for a period of time
            dataset["name"] += f", {period}-year period"
            dataset["comment"] += (
                f" Average fuel mix over a {period}"
                f"-year period {self.year}-{self.year + period}."
            )
            for exc in ws.production(dataset):
                exc["name"] += f", {period}-year period"

        # Remove existing fuel providers
        dataset["exchanges"] = [
            exc
            for exc in dataset["exchanges"]
            if exc["type"] != "technosphere"
            or (
                exc["product"] != dataset["reference product"]
                and not any(
                    x in exc["name"] for x in ["production", "evaporation", "import"]
                )
            )
        ]

        string = ""

        # if the sum is zero, we need to select a provider
        if (
            self.iam_fuel_markets.sel(region=region, variables=prod_vars)
            .interp(year=self.year)
            .sum(dim=["variables"])
            == 0
        ):
            if "hydrogen" in dataset["name"].lower():
                prod_vars = [
                    "hydrogen, from natural gas",
                ]

        for prod_var in prod_vars:
            if len(prod_vars) > 1:
                share = fuel_providers[prod_var]["find_share"](
                    prod_var, tuple(vars_map[fuel_category]), region, period
                )

            else:
                share = 1.0

            if np.isnan(share) or share <= 0:
                continue

            if isinstance(share, np.ndarray):
                share = share.item(0)

            blacklist = [
                "petroleum coke",
                "petroleum gas",
                "wax",
                "low pressure",
                "pressure, vehicle grade",
                "burned",
                "market",
            ]

            if "natural gas" in dataset["name"]:
                blacklist.remove("market")
                blacklist.append("market for natural gas, high pressure")
                blacklist.append("market group for natural gas, high pressure")

            if "low-sulfur" in dataset["name"]:
                blacklist.append("unleaded")

            if "petroleum gas" in dataset["name"]:
                blacklist.remove("petroleum gas")

            possible_names = tuple(fuel_providers[prod_var]["fuel filters"])

            possible_suppliers = self.select_multiple_suppliers(
                possible_names=possible_names,
                dataset_location=dataset["location"],
                look_for=tuple(vars_map[fuel_category]),
                blacklist=tuple(blacklist),
            )

            if not possible_suppliers:
                print(
                    f"No suppliers found for {prod_var} "
                    f"in {region} for dataset "
                    f"in location {dataset['location']}"
                )

            for supplier_key, supplier_val in possible_suppliers.items():
                # Convert m3 to kg
                conversion_factor = 1.0
                if supplier_key[-1] == "kilogram":
                    if dataset["unit"] == "cubic meter":
                        conversion_factor = 0.735

                supplier_share = share * supplier_val

                # Calculate amount of fuel input
                # Corrected by the LHV of the initial fuel
                # so that the overall composition maintains
                # the same average LHV
                amount = (
                    supplier_share
                    * (activity["lhv"] / self.fuels_specs[prod_var]["lhv"])
                    * conversion_factor
                )

                lhv = self.fuels_specs[prod_var]["lhv"] * (
                    0.735 if supplier_key[-1] == "cubic meter" else 1
                )
                co2_factor = self.fuels_specs[prod_var]["co2"] * (
                    0.735 if supplier_key[-1] == "cubic meter" else 1
                )
                biogenic_co2_share = self.fuels_specs[prod_var]["biogenic_share"]

                f_co2, nf_co2, weighted_lhv = calculate_fuel_properties(
                    amount, lhv, co2_factor, biogenic_co2_share
                )

                final_lhv += weighted_lhv
                fossil_co2 += f_co2
                non_fossil_co2 += nf_co2

                dataset = update_dataset(dataset, supplier_key, amount)

                text = (
                    f"{prod_var.capitalize()}: {(share * 100):.1f} pct @ "
                    f"{self.fuels_specs[prod_var]['lhv']} MJ/kg. "
                )
                if text not in string:
                    string += text

        if "log parameters" not in dataset:
            dataset["log parameters"] = {}

        dataset["log parameters"]["fossil CO2 per kg fuel"] = fossil_co2

        dataset["log parameters"]["non-fossil CO2 per kg fuel"] = non_fossil_co2

        dataset["log parameters"]["lower heating value"] = final_lhv

        string += f"Final average LHV of {final_lhv} MJ/kg."

        if "comment" in dataset:
            dataset["comment"] += string
        else:
            dataset["comment"] = string

        # add two new fields: `fossil CO2` and `biogenic CO2`
        dataset["fossil CO2"] = fossil_co2
        dataset["non-fossil CO2"] = non_fossil_co2
        dataset["LHV"] = final_lhv

        return dataset

    def generate_fuel_markets(self):
        """
        Create new fuel supply chains
        and update existing fuel markets.

        """

        # Create new fuel supply chains
        self.generate_fuel_supply_chains()

        fuel_markets = fetch_mapping(FUEL_MARKETS)

        # refresh the fuel filters
        # as some have been created in the meanwhile
        mapping = InventorySet(self.database)
        self.fuel_map = mapping.generate_fuel_map()
        # reverse fuel map
        self.rev_fuel_map = {}
        for fuel, activities in self.fuel_map.items():
            for activity in activities:
                self.rev_fuel_map[activity] = fuel

        d_fuels = self.get_fuel_mapping()

        vars_map = {
            "petrol, low-sulfur": [
                "petrol",
                "ethanol",
                "methanol",
                "gasoline",
                "bioethanol",
            ],
            "diesel, low-sulfur": [
                "diesel",
                "biodiesel",
            ],
            "natural gas": [
                "natural gas",
                "biomethane",
            ],
            "hydrogen": [
                "hydrogen",
            ],
            "kerosene": [
                "kerosene",
            ],
            "liquefied petroleum gas": [
                "liquefied petroleum gas",
            ],
        }

        new_datasets = []

        for fuel, activity in fuel_markets.items():
            if [
                i
                for e in self.iam_fuel_markets.variables.values
                for i in vars_map[fuel]
                if e.startswith(i)
            ]:
                prod_vars = [
                    v
                    for v in self.iam_fuel_markets.variables.values
                    if any(v.lower().startswith(i.lower()) for i in vars_map[fuel])
                ]

                d_act = self.fetch_proxies(
                    name=activity["name"],
                    ref_prod=activity["reference product"],
                    production_variable=prod_vars,
                    exact_name_match=False,
                )

                if self.system_model == "consequential":
                    periods = [
                        0,
                    ]
                else:
                    periods = [0, 20, 40, 60]

                for period in periods:
                    for region, dataset in copy.deepcopy(d_act).items():
                        for exc in ws.production(dataset):
                            if "input" in exc:
                                del exc["input"]
                        if "input" in dataset:
                            del dataset["input"]
                        if "code" in dataset:
                            dataset["code"] = str(uuid.uuid4().hex)

                        if region != "World":
                            dataset = self.generate_regional_fuel_market(
                                dataset=dataset,
                                fuel_providers=d_fuels,
                                prod_vars=prod_vars,
                                vars_map=vars_map,
                                fuel_category=fuel,
                                region=region,
                                activity=activity,
                                period=period,
                            )
                        else:
                            # World dataset
                            dataset = self.generate_world_fuel_market(
                                dataset=dataset,
                                d_act=d_act,
                                prod_vars=prod_vars,
                                period=period,
                            )

                        # add fuel market to the dictionary
                        if "log parameters" in dataset:
                            self.new_fuel_markets.update(
                                {
                                    (dataset["name"], dataset["location"]): {
                                        "fossil CO2": dataset["log parameters"][
                                            "fossil CO2 per kg fuel"
                                        ],
                                        "non-fossil CO2": dataset["log parameters"][
                                            "non-fossil CO2 per kg fuel"
                                        ],
                                        "LHV": dataset["log parameters"][
                                            "lower heating value"
                                        ],
                                    }
                                }
                            )

                        # add to log
                        self.write_log(dataset)
                        # add it to list of created datasets
                        self.add_to_index(dataset)

                        new_datasets.append(dataset)

                        # if dataset contains "low-sulfur", we make a deep copy
                        # and rename without "low-sulfur"

                        if (
                            "low-sulfur" in dataset["name"]
                            and "period" not in dataset["name"]
                        ):
                            new_dataset = copy.deepcopy(dataset)
                            new_dataset["name"] = (
                                new_dataset["name"].replace(", low-sulfur", "").strip()
                            )
                            new_dataset["reference product"] = (
                                new_dataset["reference product"]
                                .replace(", low-sulfur", "")
                                .strip()
                            )
                            for exc in ws.production(new_dataset):
                                exc["name"] = exc["name"].replace(", low-sulfur", "")
                                exc["product"] = exc["product"].replace(
                                    ", low-sulfur", ""
                                )
                            new_dataset["code"] = str(uuid.uuid4().hex)

                            # add to log
                            self.write_log(new_dataset)
                            # add it to list of created datasets
                            self.add_to_index(new_dataset)

                            new_datasets.append(new_dataset)

        # add to database
        self.database.extend(new_datasets)

        # list `market group for ` as "emptied"
        datasets_to_empty = {
            "market group for diesel",
            "market group for diesel, low-sulfur",
            "market for petrol, unleaded",
            "market for diesel",
            "market for natural gas, high pressure",
            "market group for natural gas, high pressure",
            "market for kerosene",
            "market for liquefied petroleum gas",
        }

        for old_ds in datasets_to_empty:
            for ds in ws.get_many(
                self.database,
                ws.equals("name", old_ds),
                ws.doesnt_contain_any("location", self.regions),
            ):
                self.remove_from_index(ds)

        # generate supply datasets for hydrogen
        hydrogen_supply = self.fetch_proxies(
            name="hydrogen supply, distributed by pipeline",
            ref_prod="hydrogen, gaseous, from pipeline",
            exact_name_match=True,
        )

        for region, dataset in hydrogen_supply.items():
            # replace the input of hydrogen by an input of hydrogen market
            for exc in ws.technosphere(dataset):
                if exc["name"] == "market for hydrogen, gaseous":
                    exc["location"] = region

            # add to log
            self.write_log(dataset)
            # add it to list of created datasets
            self.add_to_index(dataset)
            # add it to database
            self.database.append(dataset)

        self.relink_activities_to_new_markets()

    def write_log(self, dataset, status="created"):
        """
        Write log file.
        """

        logger.info(
            f"{status}|{self.model}|{self.scenario}|{self.year}|"
            f"{dataset['name']}|{dataset['location']}|"
            f"{dataset.get('log parameters', {}).get('initial amount of fossil CO2', '')}|"
            f"{dataset.get('log parameters', {}).get('new amount of fossil CO2', '')}|"
            f"{dataset.get('log parameters', {}).get('new amount of biogenic CO2', '')}|"
            f"{dataset.get('log parameters', {}).get('initial energy input for hydrogen production', '')}|"
            f"{dataset.get('log parameters', {}).get('new energy input for hydrogen production', '')}|"
            f"{dataset.get('log parameters', {}).get('hydrogen distribution losses', '')}|"
            f"{dataset.get('log parameters', {}).get('electricity for hydrogen compression', '')}|"
            f"{dataset.get('log parameters', {}).get('electricity for hydrogen compression after dehydrogenation', '')}|"
            f"{dataset.get('log parameters', {}).get('electricity for hydrogen pre-cooling', '')}|"
            f"{dataset.get('log parameters', {}).get('fuel conversion efficiency', '')}|"
            f"{dataset.get('log parameters', {}).get('land footprint', '')}|"
            f"{dataset.get('log parameters', {}).get('land use CO2', '')}|"
            f"{dataset.get('log parameters', {}).get('fossil CO2 per kg fuel', '')}|"
            f"{dataset.get('log parameters', {}).get('non-fossil CO2 per kg fuel', '')}|"
            f"{dataset.get('log parameters', {}).get('lower heating value', '')}"
        )
