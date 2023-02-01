"""
Integrates projections regarding fuel production and supply.
"""

import copy
import csv
from datetime import date
from pathlib import Path

import yaml
from numpy import ndarray

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
    relink_technosphere_exchanges,
    uuid,
    ws,
    wurst,
)
from .utils import DATA_DIR, get_crops_properties

REGION_CLIMATE_MAP = DATA_DIR / "fuels" / "region_to_climate.yml"
FUEL_LABELS = DATA_DIR / "fuels" / "fuel_labels.csv"
SUPPLY_CHAIN_SCENARIOS = DATA_DIR / "fuels" / "supply_chain_scenarios.yml"
HEAT_SOURCES = DATA_DIR / "fuels" / "heat_sources_map.yml"
HYDROGEN_SOURCES = DATA_DIR / "fuels" / "hydrogen_activities.yml"
HYDROGEN_SUPPLY_LOSSES = DATA_DIR / "fuels" / "hydrogen_supply_losses.yml"
METHANE_SOURCES = DATA_DIR / "fuels" / "methane_activities.yml"
LIQUID_FUEL_SOURCES = DATA_DIR / "fuels" / "liquid_fuel_activities.yml"
FUEL_MARKETS = DATA_DIR / "fuels" / "fuel_markets.yml"
BIOFUEL_SOURCES = DATA_DIR / "fuels" / "biofuels_activities.yml"


def filter_results(
    item_to_look_for: str, results: List[dict], field_to_look_at: str
) -> List[dict]:
    """Filters a list of results by a given field"""
    return [r for r in results if item_to_look_for not in r[field_to_look_at]]


def fetch_mapping(filepath: str) -> dict:
    """Returns a dictionary from a YML file"""

    with open(filepath, "r", encoding="utf-8") as stream:
        mapping = yaml.safe_load(stream)
    return mapping


def get_compression_effort(p_in: int, p_out: int, flow_rate: int) -> float:
    """Calculate the required electricity consumption from the compressor given
    an inlet and outlet pressure and a flow rate for hydrogen."""

    # result is shaft power [kW] and compressor size [kW]
    # flow_rate = mass flow rate (kg/day)
    # p_in =  input pressure (bar)
    # p_out =  output pressure (bar)
    Z_factor = 1.03198  # the hydrogen compressibility factor
    N_stages = 2  # the number of compressor stages (assumed to be 2 for this work)
    t_inlet = 310.95  # K the inlet temperature of the compressor
    y_ratio = 1.4  # the ratio of specific heats
    M_h2 = 2.15  # g/mol the molecular mass of hydrogen
    eff_comp = 0.75  # %
    R_constant = 8.314  # J/(mol*K)
    part_1 = (
        (flow_rate * (1 / (24 * 3600)))
        * ((Z_factor * t_inlet * R_constant) / (M_h2 * eff_comp))
        * ((N_stages * y_ratio / (y_ratio - 1)))
    )
    part_2 = ((p_out / p_in) ** ((y_ratio - 1) / (N_stages * y_ratio))) - 1
    power_req = part_1 * part_2
    return power_req * 24 / flow_rate


def get_pre_cooling_energy(t_amb: int, cap_util: int) -> float:
    """Calculate the required electricity consumption
    to pre-cool the hydrogen before tank filling"""

    el_pre_cooling = (0.3 / 1.6 * np.exp(-0.018 * t_amb)) + (
        (25 * np.log(t_amb) - 21) / cap_util
    )

    return el_pre_cooling


def adjust_electrolysis_electricity_requirement(year: int) -> ndarray:
    # from 58 kWh/kg H2 in 2010, down to 44 kWh in 2050
    return np.clip(-0.3538 * (year - 2010) + 58.589, 44, None)


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
    ):
        super().__init__(database, iam_data, model, pathway, year)
        # ecoinvent version
        self.version = version
        # list of fuel types
        self.fuel_labels = self.iam_data.fuel_markets.coords["variables"].values
        # dictionary of crops with specifications
        self.crops_props = get_crops_properties()
        # list to store markets that will be created
        self.new_fuel_markets = {}
        # dictionary to store mapping results, to avoid redundant effort
        self.cached_suppliers = {}

    def generate_dac_activities(self) -> None:
        """Generate regional variants of the DAC process with varying heat sources"""

        # define heat sources
        heat_map_ds = fetch_mapping(HEAT_SOURCES)

        name = "carbon dioxide, captured from atmosphere"
        original_ds = self.fetch_proxies(
            name=name,
            ref_prod="carbon dioxide",
        )

        # delete original
        self.database = [x for x in self.database if x["name"] != name]

        # loop through types of heat source
        for heat_type, activities in heat_map_ds.items():
            new_ds = copy.deepcopy(original_ds)

            new_name = f"{name}, with {heat_type}, and grid electricity"

            for _, dataset in new_ds.items():
                dataset["name"] = new_name
                dataset["code"] = str(uuid.uuid4().hex)

                for exc in ws.production(dataset):
                    exc["name"] = new_name
                    if "input" in exc:
                        del exc["input"]

                for exc in ws.technosphere(dataset):
                    if "heat" in exc["name"]:
                        exc["name"] = activities["name"]
                        exc["product"] = activities["reference product"]
                        exc["location"] = "RoW"

                        if heat_type == "heat pump heat":
                            exc["unit"] = "kilowatt hour"
                            exc["amount"] *= 1 / (
                                2.9 * 3.6
                            )  # COP of 2.9 and MJ --> kWh
                            exc["location"] = "RER"

                            dataset[
                                "comment"
                            ] = "Dataset generated by `premise`, initially based on Terlouw et al. 2021. "
                            dataset["comment"] += (
                                "A CoP of 2.9 is assumed for the heat pump. But the heat pump itself is not"
                                + " considered here. "
                            )

                dataset["comment"] += (
                    "The CO2 is compressed from 1 bar to 25 bar, "
                    + " for which 0.78 kWh is considered. Furthermore, there's a 2.1% loss on site"
                    + " and only a 1 km long pipeline transport."
                )

                self.cache, dataset = relink_technosphere_exchanges(
                    dataset,
                    self.database,
                    self.model,
                    cache=self.cache,
                    contained=False,
                )
            self.database.extend(new_ds.values())

            # Add created dataset to `self.list_datasets`
            self.list_datasets.extend(
                [
                    (
                        act["name"],
                        act["reference product"],
                        act["location"],
                    )
                    for act in new_ds.values()
                ]
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
        self, name: str, ref_prod: str, unit: str, loc: str, exclude: List[str] = None
    ) -> Dict[Tuple[Any, Any, Any, Any], float]:
        """
        Return a list of potential suppliers given a name, reference product,
        unit and location, with their respective supply share (based on production volumes).

        :param name:
        :param ref_prod:
        :param unit:
        :param loc:
        :param exclude:
        :return:
        """

        # if we find a result in the cache dictionary, return it
        if (name, ref_prod, loc) in self.cached_suppliers:
            return self.cached_suppliers[(name, ref_prod, loc)]

        ecoinvent_regions = self.geo.iam_to_ecoinvent_location(loc)
        # search first for suppliers in `loc`, but also in the ecoinvent
        # locations to are comprised in `loc`, and finally in `RER`, `RoW` and `GLO`,
        possible_locations = [[loc], ecoinvent_regions, ["RER"], ["RoW"], ["GLO"]]
        suppliers, counter = [], 0

        # while we do nto find a result
        while len(suppliers) == 0:
            suppliers = list(
                get_suppliers_of_a_region(
                    database=self.database,
                    locations=possible_locations[counter],
                    names=[name] if isinstance(name, str) else name,
                    reference_product=ref_prod,
                    unit=unit,
                    exclude=exclude,
                )
            )
            counter += 1

        # find production volume-based share
        suppliers = get_shares_from_production_volume(suppliers)

        # store the result in cache for next time
        self.cached_suppliers[(name, ref_prod, loc)] = suppliers

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

        We also assume efficiency gains over time for the PEM electrolysis process: from 58 kWh/kg H2 in 2010,
        down to 44 kWh by 2050, according to a literature review conducted by the Paul Scherrer Institute.

        """

        hydrogen_sources = fetch_mapping(HYDROGEN_SOURCES)

        for hydrogen_type, hydrogen_activity in hydrogen_sources.items():
            new_ds = self.fetch_proxies(
                name=hydrogen_activity["name"],
                ref_prod="Hydrogen",
                production_variable=hydrogen_activity["var"],
                relink=True,
            )

            for region, dataset in new_ds.items():
                # we adjust the electrolysis efficiency
                if hydrogen_type == "from electrolysis":
                    for exc in ws.technosphere(dataset):
                        if "market group for electricity" in exc["name"]:
                            exc["amount"] = adjust_electrolysis_electricity_requirement(
                                self.year
                            )

                    string = f" The electricity input per kg of H2 has been adapted to the year {self.year}."
                    if "comment" in dataset:
                        dataset["comment"] += string
                    else:
                        dataset["comment"] = string

                dataset[
                    "comment"
                ] = "Region-specific hydrogen production dataset generated by `premise`. "

            self.database.extend(new_ds.values())

        print("Generate region-specific hydrogen supply chains.")

        # loss coefficients for hydrogen supply
        losses = fetch_mapping(HYDROGEN_SUPPLY_LOSSES)

        supply_chain_scenarios = fetch_mapping(SUPPLY_CHAIN_SCENARIOS)

        for act in [
            "hydrogen embrittlement inhibition",
            "geological hydrogen storage",
            "hydrogenation of hydrogen",
            "dehydrogenation of hydrogen",
            "Hydrogen refuelling station",
        ]:
            new_ds = self.fetch_proxies(name=act, ref_prod=" ")

            for _, dataset in new_ds.items():
                for exc in ws.production(dataset):
                    if "input" in exc:
                        del exc["input"]

                self.cache, dataset = relink_technosphere_exchanges(
                    dataset, self.database, self.model, cache=self.cache
                )

            self.database.extend(new_ds.values())

        for region in self.regions:
            for hydrogen_type, hydrogen_activity in hydrogen_sources.items():
                for vehicle, config in supply_chain_scenarios.items():
                    for state in config["state"]:
                        for distance in config["distance"]:
                            # dataset creation
                            new_act = {
                                "location": region,
                                "name": f"hydrogen supply, {hydrogen_type}, by {vehicle}, as {state}, over {distance} km",
                                "reference product": "hydrogen, 700 bar",
                                "unit": "kilogram",
                                "database": self.database[1]["database"],
                                "code": str(uuid.uuid4().hex),
                                "comment": "Dataset representing hydrogen supply, generated by `premise`.",
                            }

                            # production flow
                            new_exc = [
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
                            ]

                            # transport
                            for trspt in config["vehicle"]:
                                suppliers = self.find_suppliers(
                                    name=trspt["name"],
                                    ref_prod=trspt["reference product"],
                                    unit=trspt["unit"],
                                    loc=region,
                                )

                                for supplier, share in suppliers.items():
                                    new_exc.append(
                                        {
                                            "uncertainty type": 0,
                                            "amount": distance * share / 1000
                                            if supplier[-1] == "ton kilometer"
                                            else distance
                                            * share
                                            / 2
                                            * (1 / eval(config["lifetime"])),
                                            "type": "technosphere",
                                            "product": supplier[2],
                                            "name": supplier[0],
                                            "unit": supplier[-1],
                                            "location": supplier[1],
                                            "comment": f"Transport over {distance} km by {vehicle}. ",
                                        }
                                    )

                                string = f"Transport over {distance} km by {vehicle}."

                                if "comment" in new_act:
                                    new_act["comment"] += string
                                else:
                                    new_act["comment"] = string

                            # need for inhibitor and purification if CNG pipeline
                            # electricity for purification: 2.46 kWh/kg H2
                            if vehicle == "CNG pipeline":
                                inhibbitor_ds = ws.get_one(
                                    self.database,
                                    ws.contains(
                                        "name", "hydrogen embrittlement inhibition"
                                    ),
                                    ws.equals("location", region),
                                )

                                new_exc.append(
                                    {
                                        "uncertainty type": 0,
                                        "amount": 1,
                                        "type": "technosphere",
                                        "product": inhibbitor_ds["reference product"],
                                        "name": inhibbitor_ds["name"],
                                        "unit": inhibbitor_ds["unit"],
                                        "location": region,
                                        "comment": "Injection of an inhibiting gas (oxygen) "
                                        "to prevent embritllement of metal. ",
                                    }
                                )

                                string = (
                                    "2.46 kWh/kg H2 is needed to purify the hydrogen from the inhibiting gas. "
                                    "The recovery rate for hydrogen after separation from the inhibitor gas is 93%. "
                                )
                                if "comment" in new_act:
                                    new_act["comment"] += string
                                else:
                                    new_act["comment"] = string

                            if "regional storage" in config:
                                storage_ds = ws.get_one(
                                    self.database,
                                    ws.contains(
                                        "name", config["regional storage"]["name"]
                                    ),
                                    ws.contains(
                                        "reference product",
                                        config["regional storage"]["reference product"],
                                    ),
                                    ws.equals("location", region),
                                    ws.equals(
                                        "unit", config["regional storage"]["unit"]
                                    ),
                                )

                                new_exc.append(
                                    {
                                        "uncertainty type": 0,
                                        "amount": 1,
                                        "type": "technosphere",
                                        "product": storage_ds["reference product"],
                                        "name": storage_ds["name"],
                                        "unit": storage_ds["unit"],
                                        "location": region,
                                        "comment": "Geological storage (salt cavern).",
                                    }
                                )

                                string = (
                                    " Geological storage is added. It includes 0.344 kWh for "
                                    "the injection and pumping of 1 kg of H2."
                                )
                                if "comment" in new_act:
                                    new_act["comment"] += string
                                else:
                                    new_act["comment"] = string

                            # electricity for compression
                            if state in ["gaseous", "liquid"]:
                                # if gaseous
                                # if transport by truck, compression from 25 bar to 500 bar for the transport
                                # and from 500 bar to 900 bar for dispensing in 700 bar storage tanks

                                # if transport by pipeline, initial compression from 25 bar to 100 bar
                                # and 0.6 kWh re-compression every 250 km
                                # and finally from 100 bar to 900 bar for dispensing in 700 bar storage tanks

                                # if liquid
                                # liquefaction electricity need
                                # currently, 12 kWh/kg H2
                                # mid-term, 8 kWh/ kg H2
                                # by 2050, 6 kWh/kg H2

                                if state == "gaseous":
                                    if vehicle == "truck":
                                        electricity_comp = get_compression_effort(
                                            25, 500, 1000
                                        )
                                        electricity_comp += get_compression_effort(
                                            500, 900, 1000
                                        )
                                    else:
                                        electricity_comp = get_compression_effort(
                                            25, 100, 1000
                                        ) + (0.6 * distance / 250)
                                        electricity_comp += get_compression_effort(
                                            100, 900, 1000
                                        )

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

                                if "comment" in new_act:
                                    new_act["comment"] += string
                                else:
                                    new_act["comment"] = string

                            # electricity for hydrogenation, dehydrogenation and
                            # compression at delivery
                            if state == "liquid organic compound":
                                hydrogenation_ds = ws.get_one(
                                    self.database,
                                    ws.equals("name", "hydrogenation of hydrogen"),
                                    ws.equals("location", region),
                                )

                                dehydrogenation_ds = ws.get_one(
                                    self.database,
                                    ws.equals("name", "dehydrogenation of hydrogen"),
                                    ws.equals("location", region),
                                )

                                new_exc.extend(
                                    [
                                        {
                                            "uncertainty type": 0,
                                            "amount": 1,
                                            "type": "technosphere",
                                            "product": hydrogenation_ds[
                                                "reference product"
                                            ],
                                            "name": hydrogenation_ds["name"],
                                            "unit": hydrogenation_ds["unit"],
                                            "location": region,
                                        },
                                        {
                                            "uncertainty type": 0,
                                            "amount": 1,
                                            "type": "technosphere",
                                            "product": dehydrogenation_ds[
                                                "reference product"
                                            ],
                                            "name": dehydrogenation_ds["name"],
                                            "unit": dehydrogenation_ds["unit"],
                                            "location": region,
                                        },
                                    ]
                                )

                                # After dehydrogenation at ambient temperature at delivery
                                # the hydrogen needs to be compressed up to 900 bar to be dispensed
                                # in 700 bar storage tanks

                                electricity_comp = get_compression_effort(25, 900, 1000)

                                suppliers = self.find_suppliers(
                                    name="market group for electricity, low voltage",
                                    ref_prod="electricity, low voltage",
                                    unit="kilowatt hour",
                                    loc=region,
                                    exclude=["period"],
                                )

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

                                string = (
                                    " Hydrogenation and dehydrogenation of hydrogen included. "
                                    "Compression at delivery after dehydrogenation also included."
                                )
                                if "comment" in new_act:
                                    new_act["comment"] += string
                                else:
                                    new_act["comment"] = string

                            # fetch the H2 production activity
                            h2_ds = ws.get_one(
                                self.database,
                                ws.equals("name", hydrogen_activity["name"]),
                                ws.equals("location", region),
                            )

                            # include losses along the way
                            string = ""
                            total_loss = 1
                            for loss, val in losses[vehicle][state].items():
                                val = float(val)

                                if loss == "boil-off":
                                    if vehicle == "truck":
                                        # average truck speed
                                        speed = 50
                                    else:
                                        # average ship speed
                                        speed = 36
                                    days = distance / speed / 24
                                    # boil-off losses, function of days in transit
                                    total_loss *= np.power(1 + val, days)

                                    string += f"Boil-off losses: {int((np.power(1 + val, days) - 1) * 100)}%. "

                                elif loss == "pipeline_leak":
                                    # pipeline losses, function of distance
                                    total_loss *= 1 + (val * distance)
                                    string += f"Pipeline losses: {int((1 + (val * distance) - 1) * 100)}%. "
                                else:
                                    total_loss *= 1 + val
                                    string += f"{loss} losses: {int(val * 100)}%. "

                            new_exc.append(
                                {
                                    "uncertainty type": 0,
                                    "amount": 1 * total_loss,
                                    "type": "technosphere",
                                    "product": h2_ds["reference product"],
                                    "name": h2_ds["name"],
                                    "unit": h2_ds["unit"],
                                    "location": region,
                                }
                            )

                            # adds losses as hydrogen emission to air
                            new_exc.append(
                                {
                                    "uncertainty type": 0,
                                    "amount": total_loss - 1,
                                    "type": "biosphere",
                                    "name": "Hydrogen",
                                    "unit": "kilogram",
                                    "categories": ("air",),
                                    "input": (
                                        "biosphere3",
                                        "b301fa9a-ba60-4eac-8ccc-6ccbdf099b35",
                                    ),
                                }
                            )

                            if "comment" in new_act:
                                new_act["comment"] += string
                            else:
                                new_act["comment"] = string

                            # add fuelling station, including storage tank
                            ds_h2_station = ws.get_one(
                                self.database,
                                ws.equals("name", "Hydrogen refuelling station"),
                                ws.equals("location", region),
                            )

                            new_exc.append(
                                {
                                    "uncertainty type": 0,
                                    "amount": 1
                                    / (
                                        600 * 365 * 40
                                    ),  # 1 over lifetime: 40 years, 600 kg H2/day
                                    "type": "technosphere",
                                    "product": ds_h2_station["reference product"],
                                    "name": ds_h2_station["name"],
                                    "unit": ds_h2_station["unit"],
                                    "location": region,
                                }
                            )

                            # finally, add pre-cooling
                            # pre-cooling is needed before filling vehicle tanks
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
                            cap_util = np.interp(self.year, [2020, 2050], [10, 150])
                            el_pre_cooling = get_pre_cooling_energy(t_amb, cap_util)

                            suppliers = self.find_suppliers(
                                name="market group for electricity, low voltage",
                                ref_prod="electricity, low voltage",
                                unit="kilowatt hour",
                                loc=region,
                                exclude=["period"],
                            )

                            for supplier, share in suppliers.items():
                                new_exc.append(
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
                            if "comment" in new_act:
                                new_act["comment"] += string
                            else:
                                new_act["comment"] = string

                            new_act["exchanges"] = new_exc

                            self.cache, new_act = relink_technosphere_exchanges(
                                new_act, self.database, self.model, cache=self.cache
                            )

                            self.database.append(new_act)

                            # Add created dataset to `self.list_datasets`
                            self.list_datasets.append(
                                (
                                    new_act["name"],
                                    new_act["reference product"],
                                    new_act["location"],
                                )
                            )

    def generate_biogas_activities(self):
        """
        Generate biogas activities.
        """

        fuel_activities = fetch_mapping(METHANE_SOURCES)

        for fuel, activities in fuel_activities.items():
            for activity in activities:
                original_ds = self.fetch_proxies(name=activity, ref_prod=" ")

                # delete original
                # self.database = [x for x in self.database if x["name"] != f]

                if fuel == "methane, synthetic":
                    for co2_type in [
                        (
                            "carbon dioxide, captured from atmosphere, with waste heat, and grid electricity",
                            "carbon dioxide, captured from the atmosphere",
                            "waste heat",
                        ),
                        (
                            "carbon dioxide, captured from atmosphere, with heat pump heat, and grid electricity",
                            "carbon dioxide, captured from the atmosphere",
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
                                    exc["name"] = co2_type[0]
                                    exc["product"] = co2_type[1]
                                    exc["location"] = region

                            for exc in ws.technosphere(dataset):
                                if (
                                    "methane, from electrochemical methanation"
                                    in exc["name"]
                                ):
                                    exc["name"] += f", using {co2_type[2]}"
                                    exc["location"] = dataset["location"]

                                    dataset["name"] = dataset["name"].replace(
                                        "from electrochemical methanation",
                                        f"from electrochemical methanation (H2 from electrolysis, CO2 from DAC using {co2_type[2]})",
                                    )

                                    for prod in ws.production(dataset):
                                        prod["name"] = prod["name"].replace(
                                            "from electrochemical methanation",
                                            f"from electrochemical methanation (H2 from electrolysis, CO2 from DAC using {co2_type[2]})",
                                        )

                        self.database.extend(new_ds.values())

                        # Add created dataset to `self.list_datasets`
                        self.list_datasets.extend(
                            [
                                (
                                    act["name"],
                                    act["reference product"],
                                    act["location"],
                                )
                                for act in new_ds.values()
                            ]
                        )

                else:
                    new_ds = copy.deepcopy(original_ds)

                    for region, dataset in new_ds.items():
                        dataset["code"] = str(uuid.uuid4().hex)
                        for exc in ws.production(dataset):
                            if "input" in exc:
                                exc.pop("input")

                        self.cache, dataset = relink_technosphere_exchanges(
                            dataset, self.database, self.model, cache=self.cache
                        )

                    self.database.extend(new_ds.values())

                    # Add created dataset to `self.list_datasets`
                    self.list_datasets.extend(
                        [
                            (
                                act["name"],
                                act["reference product"],
                                act["location"],
                            )
                            for act in new_ds.values()
                        ]
                    )

    def generate_synthetic_fuel_activities(self):
        """
        Generate synthetic fuel activities.
        """

        fuel_activities = fetch_mapping(LIQUID_FUEL_SOURCES)

        for activities in fuel_activities.values():
            for activity in activities:
                new_ds = self.fetch_proxies(name=activity, ref_prod=" ")
                for region, dataset in new_ds.items():
                    for exc in ws.production(dataset):
                        if "input" in exc:
                            del exc["input"]

                    for exc in ws.technosphere(dataset):
                        if "carbon dioxide, captured from atmosphere" in exc["name"]:
                            exc[
                                "name"
                            ] = "carbon dioxide, captured from atmosphere, with heat pump heat, and grid electricity"
                            exc[
                                "product"
                            ] = "carbon dioxide, captured from the atmosphere"
                            exc["location"] = region
                            exc["unit"] = "kilogram"

                    self.cache, dataset = relink_technosphere_exchanges(
                        dataset, self.database, self.model, cache=self.cache
                    )

                    self.database.append(dataset)

    def adjust_land_use(self, dataset, region, crop_type):
        """
        Adjust land use.
        """
        for exc in dataset["exchanges"]:
            # we adjust the land use
            if exc["type"] == "biosphere" and exc["name"].startswith("Occupation"):
                # lower heating value, dry basis
                lhv_ar = dataset["LHV [MJ/kg dry]"]

                # Ha/GJ
                land_use = (
                    self.iam_data.land_use.sel(region=region, variables=crop_type)
                    .interp(year=self.year)
                    .values
                )

                if land_use > 0:
                    # HA to m2
                    land_use *= 10000
                    # m2/GJ to m2/MJ
                    land_use /= 1000
                    # m2/kg, as received
                    land_use *= lhv_ar
                    # update exchange value
                    exc["amount"] = float(land_use)

                    string = (
                        f"The land area occupied has been modified to {land_use}, "
                        f"to be in line with the scenario {self.scenario} of {self.model.upper()} "
                        f"in {self.year} in the region {region}. "
                    )
                    if "comment" in dataset:
                        dataset["comment"] += string
                    else:
                        dataset["comment"] = string

        return dataset

    def adjust_land_use_change_emissions(self, dataset, region, crop_type):
        """
        Adjust land use change emissions.
        """

        # then, we should include the Land Use Change-induced CO2 emissions
        # those are given in kg CO2-eq./GJ of primary crop energy

        # kg CO2/GJ
        land_use_co2 = (
            self.iam_data.land_use_change.sel(region=region, variables=crop_type)
            .interp(year=self.year)
            .values
        )

        if land_use_co2 > 0:
            # lower heating value, as received
            lhv_ar = dataset["LHV [MJ/kg dry]"]

            # kg CO2/MJ
            land_use_co2 /= 1000
            land_use_co2 *= lhv_ar

            land_use_co2_exc = {
                "uncertainty type": 0,
                "loc": float(land_use_co2),
                "amount": float(land_use_co2),
                "type": "biosphere",
                "name": "Carbon dioxide, from soil or biomass stock",
                "unit": "kilogram",
                "input": (
                    "biosphere3",
                    "78eb1859-abd9-44c6-9ce3-f3b5b33d619c",
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

        return dataset

    def adjust_biomass_conversion_efficiency(self, dataset, region, crop_type):
        """
        Adjust biomass conversion efficiency.
        """

        # modify efficiency
        # fetch the `scaling factor`, compared to 2020
        variables = [
            v
            for v in self.iam_data.efficiency.variables.values
            if crop_type.lower() in v.lower()
            and any(x.lower() in v.lower() for x in ["bioethanol", "biodiesel"])
        ]

        if variables:
            scaling_factor = 1 / self.find_iam_efficiency_change(
                variable=variables, location=region
            )

            # Rescale the biomass input according to the IAM efficiency values
            wurst.change_exchanges_by_constant_factor(
                dataset,
                scaling_factor,
                [
                    ws.equals("unit", "kilogram"),
                    ws.either(
                        ws.contains("name", "Farming"), ws.contains("name", "Supply")
                    ),
                ],
                [ws.doesnt_contain_any("name", self.emissions_map)],
            )

            string = (
                f"The biomass input has been rescaled by `premise` "
                f"by {int((scaling_factor - 1) * 100)}%.\n"
                f"To be in line with the scenario {self.scenario} of {self.model.upper()} "
                f"in {self.year} in the region {region}.\n"
            )

            if "comment" in dataset:
                dataset["comment"] += string
            else:
                dataset["comment"] = string

            if "ethanol" in dataset["name"].lower():
                dataset[
                    "comment"
                ] += (
                    "Bioethanol has a combustion CO2 emission factor of 1.91 kg CO2/kg."
                )
            if "biodiesel" in dataset["name"].lower():
                dataset[
                    "comment"
                ] += "Biodiesel has a combustion CO2 emission factor of 2.85 kg CO2/kg."

        return dataset

    def generate_biofuel_activities(self):
        """
        Create region-specific biofuel datasets.
        Update the conversion efficiency.
        :return:
        """

        # region -> climate dictionary
        d_region_climate = fetch_mapping(REGION_CLIMATE_MAP)[self.model]
        # climate --> {crop type --> crop} dictionary
        crop_types = list(self.crops_props.keys())
        climates = set(d_region_climate.values())
        d_climate_crop_type = {
            clim: {
                crop_type: self.crops_props[crop_type]["crop_type"][self.model][clim]
                for crop_type in crop_types
            }
            for clim in climates
        }

        biofuel_activities = fetch_mapping(BIOFUEL_SOURCES)
        l_crops = []

        for climate in ["tropical", "temperate"]:
            regions = [k for k, v in d_region_climate.items() if v == climate]
            for crop_type in d_climate_crop_type[climate]:
                specific_crop = d_climate_crop_type[climate][crop_type]
                names = biofuel_activities[crop_type][specific_crop]

                if specific_crop not in l_crops:
                    l_crops.append(specific_crop)
                else:
                    continue

                if specific_crop == "corn":
                    regions = list(d_region_climate.keys())

                for name in names:
                    try:
                        prod_label = [
                            l
                            for l in self.fuel_labels
                            if crop_type.lower() in l.lower()
                            and any(
                                i.lower() in l.lower()
                                for i in ("biodiesel", "bioethanol")
                            )
                        ][0]
                    except IndexError:
                        continue

                    new_ds = self.fetch_proxies(
                        name=name,
                        ref_prod=" ",
                        production_variable=prod_label,
                        regions=regions,
                    )

                    for region in regions:
                        # if this is a fuel production activity
                        # we need to adjust the process efficiency
                        if any(
                            i in new_ds[region]["name"]
                            for i in ["Ethanol production", "Biodiesel production"]
                        ):
                            new_ds[region] = self.adjust_biomass_conversion_efficiency(
                                dataset=new_ds[region],
                                region=region,
                                crop_type=crop_type,
                            )

                        # if this is a farming activity
                        # and if the product (crop) is not a residue
                        # and if we have land use info from the IAM
                        if self.iam_data.land_use is not None:
                            if (
                                "farming and supply" in new_ds[region]["name"].lower()
                                and crop_type.lower()
                                in self.iam_data.land_use.variables.values
                                and not any(
                                    i in new_ds[region]["name"].lower()
                                    for i in ["straw", "residue", "stover"]
                                )
                            ):
                                new_ds[region] = self.adjust_land_use(
                                    dataset=new_ds[region],
                                    region=region,
                                    crop_type=crop_type,
                                )

                        # if this is a farming activity
                        # and if the product (crop) is not a residue
                        # and if we have land use change CO2 info from the IAM
                        if not self.iam_data.land_use_change is None:
                            if (
                                "farming and supply" in new_ds[region]["name"].lower()
                                and crop_type.lower()
                                in self.iam_data.land_use_change.variables.values
                                and not any(
                                    i in new_ds[region]["name"].lower()
                                    for i in ["straw", "residue", "stover"]
                                )
                            ):
                                new_ds[region] = self.adjust_land_use_change_emissions(
                                    dataset=new_ds[region],
                                    region=region,
                                    crop_type=crop_type,
                                )

                    self.database.extend(new_ds.values())

    def get_fuel_mapping(self):
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
            for fuel in self.iam_data.fuel_markets.variables.values
        }

    def fetch_fuel_share(self, fuel, fuel_types, region):
        """Return a fuel mix for a specific IAM region, for a specific year."""

        variables = [
            v
            for v in self.iam_data.fuel_markets.variables.values
            if any(x.lower() in v.lower() for x in fuel_types)
        ]

        return (
            self.iam_data.fuel_markets.sel(region=region, variables=fuel)
            / self.iam_data.fuel_markets.sel(region=region, variables=variables).sum(
                dim="variables"
            )
        ).values

    def relink_activities_to_new_markets(self):
        """
        Links fuel input exchanges to new datasets with the appropriate IAM location.

        Does not return anything.
        """

        # Filter all activities that consume fuels
        acts_to_ignore = set(x[0] for x in self.new_fuel_markets)

        fuel_markets = fetch_mapping(FUEL_MARKETS)

        for dataset in ws.get_many(
            self.database,
            ws.exclude(ws.either(*[ws.contains("name", n) for n in acts_to_ignore])),
        ):
            # check that a fuel input exchange is present in the list of inputs
            # check also for "market group for" inputs
            if any(
                f[0] == exc["name"]
                for exc in dataset["exchanges"]
                for f in self.new_fuel_markets
            ) or any(
                exc["name"] == f[0].replace("market for", "market group for")
                for exc in dataset["exchanges"]
                for f in self.new_fuel_markets
            ):
                amount_fossil_co2, amount_non_fossil_co2 = [0, 0]

                for _, activity in fuel_markets.items():
                    # checking that it is one of the markets
                    # that has been newly created
                    if activity["name"] in acts_to_ignore:
                        excs = list(
                            ws.get_many(
                                dataset["exchanges"],
                                ws.either(
                                    ws.equals("name", activity["name"]),
                                    ws.equals(
                                        "name",
                                        activity["name"].replace(
                                            "market for", "market group for"
                                        ),
                                    ),
                                ),
                                ws.equals("unit", activity["unit"]),
                                ws.equals("type", "technosphere"),
                            )
                        )

                        amount = 0
                        for exc in excs:
                            amount += exc["amount"]
                            dataset["exchanges"].remove(exc)

                        if amount > 0:
                            if dataset["location"] in self.regions:
                                supplier_loc = dataset["location"]

                            else:
                                new_loc = self.geo.ecoinvent_to_iam_location(
                                    dataset["location"]
                                )
                                supplier_loc = (
                                    new_loc
                                    if new_loc in self.regions
                                    else self.regions[0]
                                )

                            new_exc = {
                                "name": activity["name"],
                                "product": activity["reference product"],
                                "amount": amount,
                                "type": "technosphere",
                                "unit": activity["unit"],
                                "location": supplier_loc,
                            }

                            dataset["exchanges"].append(new_exc)

                            amount_fossil_co2 += (
                                amount
                                * self.new_fuel_markets[
                                    (activity["name"], supplier_loc)
                                ]["fossil CO2"]
                            )
                            amount_non_fossil_co2 += (
                                amount
                                * self.new_fuel_markets[
                                    (activity["name"], supplier_loc)
                                ]["non-fossil CO2"]
                            )

                # update fossil and biogenic CO2 emissions
                list_items_to_ignore = [
                    "blending",
                    "market group",
                    "lubricating oil production",
                    "petrol production",
                ]
                if amount_non_fossil_co2 > 0 and not any(
                    x in dataset["name"].lower() for x in list_items_to_ignore
                ):
                    # test for the presence of a fossil CO2 flow
                    if not [
                        e
                        for e in dataset["exchanges"]
                        if "Carbon dioxide, fossil" in e["name"]
                    ]:
                        print(
                            f"Warning: {dataset['name'], dataset['location']} has no fossil CO2 output flow."
                        )

                    # subtract the biogenic CO2 amount to the
                    # initial fossil CO2 emission amount

                    for exc in ws.biosphere(
                        dataset, ws.equals("name", "Carbon dioxide, fossil")
                    ):
                        if (exc["amount"] - amount_non_fossil_co2) < 0:
                            exc["amount"] = 0
                        else:
                            exc["amount"] -= amount_non_fossil_co2

                    # add the biogenic CO2 emission flow
                    non_fossil_co2 = {
                        "uncertainty type": 0,
                        "loc": amount_non_fossil_co2,
                        "amount": amount_non_fossil_co2,
                        "type": "biosphere",
                        "name": "Carbon dioxide, non-fossil",
                        "unit": "kilogram",
                        "categories": ("air",),
                        "input": ("biosphere3", "eba59fd6-f37e-41dc-9ca3-c7ea22d602c7"),
                    }

                    dataset["exchanges"].append(non_fossil_co2)

    def select_multiple_suppliers(self, fuel, d_fuels, dataset, look_for):
        """
        Select multiple suppliers for a specific fuel.
        """

        # We have several potential fuel suppliers
        # We will look up their respective production volumes
        # And include them proportionally to it

        ecoinvent_regions = self.geo.iam_to_ecoinvent_location(dataset["location"])
        possible_locations = [
            dataset["location"],
            [*ecoinvent_regions],
            "RER",
            "Europe without Switzerland",
            "RoW",
            "GLO",
        ]
        possible_names = list(d_fuels[fuel]["fuel filters"])
        suppliers, counter = [], 0

        while not suppliers:
            suppliers = list(
                ws.get_many(
                    self.database,
                    ws.either(*[ws.contains("name", sup) for sup in possible_names]),
                    ws.either(
                        *[
                            ws.equals("location", item)
                            for item in possible_locations[counter]
                        ]
                    )
                    if isinstance(possible_locations[counter], list)
                    else ws.equals("location", possible_locations[counter]),
                    ws.either(
                        *[ws.contains("reference product", item) for item in look_for]
                    ),
                    ws.doesnt_contain_any(
                        "reference product",
                        [
                            "petroleum coke",
                            "petroleum gas",
                            "wax",
                            "low pressure",
                            "pressure, vehicle grade",
                        ],
                    ),
                )
            )
            counter += 1

            if "low-sulfur" in dataset["name"]:
                suppliers = filter_results(
                    item_to_look_for="unleaded",
                    results=suppliers,
                    field_to_look_at="reference product",
                )

            if "unleaded" in dataset["name"]:
                suppliers = filter_results(
                    item_to_look_for="low-sulfur",
                    results=suppliers,
                    field_to_look_at="reference product",
                )

        suppliers = get_shares_from_production_volume(suppliers)

        return suppliers

    def generate_fuel_supply_chains(self):
        """Duplicate fuel chains and make them IAM region-specific"""

        # DAC datasets
        print("Generate region-specific direct air capture processes.")
        self.generate_dac_activities()

        # hydrogen
        print("Generate region-specific hydrogen production pathways.")
        self.generate_hydrogen_activities()

        # biogas
        print("Generate region-specific biogas and syngas supply chains.")
        self.generate_biogas_activities()

        # synthetic fuels
        print("Generate region-specific synthetic fuel supply chains.")
        self.generate_synthetic_fuel_activities()

        # biofuels
        print("Generate region-specific biofuel supply chains.")
        self.generate_biofuel_activities()

    def generate_fuel_markets(self):
        """Create new fuel supply chains
        and update existing fuel markets"""

        # Create new fuel supply chains
        self.generate_fuel_supply_chains()

        print("Generate new fuel markets.")

        # we start by creating region-specific "diesel, burned in" markets
        prod_vars = [p for p in self.fuel_labels if "diesel" in p]
        datasets_to_create = [
            (
                "diesel, burned in agricultural machinery",
                "diesel, burned in agricultural machinery",
            ),
            (
                "diesel, burned in building machine",
                "diesel, burned in building machine",
            ),
            (
                "diesel, burned in diesel-electric generating set",
                "diesel, burned in diesel-electric generating set",
            ),
            (
                "diesel, burned in diesel-electric generating set, 10MW",
                "diesel, burned in diesel-electric generating set, 10MW",
            ),
            (
                "diesel, burned in diesel-electric generating set, 18.5kW",
                "diesel, burned in diesel-electric generating set, 18.5kW",
            ),
            ("diesel, burned in fishing vessel", "diesel, burned in fishing vessel"),
            (
                "market for diesel, burned in agricultural machinery",
                "diesel, burned in agricultural machinery",
            ),
            (
                "market for diesel, burned in building machine",
                "diesel, burned in building machine",
            ),
            (
                "market for diesel, burned in diesel-electric generating set, 10MW",
                "diesel, burned in diesel-electric generating set, 10MW",
            ),
            (
                "market for diesel, burned in diesel-electric generating set, 18.5kW",
                "diesel, burned in diesel-electric generating set, 18.5kW",
            ),
            (
                "market for diesel, burned in fishing vessel",
                "diesel, burned in fishing vessel",
            ),
        ]

        for dataset in datasets_to_create:
            new_ds = self.fetch_proxies(
                name=dataset[0],
                ref_prod=dataset[1],
                production_variable=prod_vars,
                relink=True,
            )

            self.database.extend(new_ds.values())

        fuel_markets = fetch_mapping(FUEL_MARKETS)

        # refresh the fuel filters
        # as some have been created in the meanwhile
        mapping = InventorySet(self.database)
        self.fuel_map = mapping.generate_fuel_map()
        d_fuels = self.get_fuel_mapping()

        # to log new fuel markets
        new_fuel_markets = [
            [
                "market name",
                "location",
                "unit",
                "reference product",
                "fuel type",
                "supplier name",
                "supplier reference product",
                "supplier location",
                "supplier unit",
                "fuel mix share (energy-wise)",
                "amount supplied [kg]",
                "LHV [mj/kg]",
                "CO2 emmission factor [kg CO2]",
                "biogenic share",
            ]
        ]

        vars_map = {
            "petrol, unleaded": ["petrol", "ethanol", "methanol", "gasoline"],
            "petrol, low-sulfur": ["petrol", "ethanol", "methanol", "gasoline"],
            "diesel, low-sulfur": ["diesel", "biodiesel"],
            "diesel": ["diesel", "biodiesel"],
            "natural gas": ["natural gas", "biomethane"],
            "hydrogen": ["hydrogen"],
        }

        for fuel, activity in fuel_markets.items():
            if [i for e in self.fuel_labels for i in vars_map[fuel] if i in e]:
                print(f"--> {fuel}")

                prod_vars = [
                    v
                    for v in self.iam_data.fuel_markets.variables.values
                    if any(i.lower() in v.lower() for i in vars_map[fuel])
                ]

                d_act = self.fetch_proxies(
                    name=activity["name"],
                    ref_prod=activity["reference product"],
                    production_variable=prod_vars,
                )

                for region, dataset in d_act.items():
                    if region != "World":
                        string = " Fuel market composition: "
                        fossil_co2, non_fossil_co2, final_lhv = [0, 0, 0]

                        # remove existing fuel providers
                        dataset["exchanges"] = [
                            exc
                            for exc in dataset["exchanges"]
                            if exc["type"] != "technosphere"
                            or (
                                exc["product"] != dataset["reference product"]
                                and not any(
                                    x in exc["name"]
                                    for x in ["production", "evaporation", "import"]
                                )
                            )
                        ]

                        for var in prod_vars:
                            share = d_fuels[var]["find_share"](
                                var, vars_map[fuel], region
                            )

                            if share > 0:
                                possible_suppliers = self.select_multiple_suppliers(
                                    var, d_fuels, dataset, vars_map[fuel]
                                )
                                if not possible_suppliers:
                                    print(
                                        f"ISSUE with {var} in {region} for ds in location {dataset['location']}"
                                    )

                                for key, val in possible_suppliers.items():
                                    # m3 to kg conversion
                                    if key[-1] != activity["unit"]:
                                        conversion_factor = 0.679
                                    else:
                                        # kg to kg
                                        conversion_factor = 1

                                    supplier_share = share * val

                                    # LHV of the fuel before update
                                    reference_lhv = activity["lhv"]

                                    # amount of fuel input
                                    # corrected by the LHV of the initial fuel
                                    # so that the overall composition maintains
                                    # the same average LHV

                                    amount = (
                                        supplier_share
                                        * (reference_lhv / self.fuels_specs[var]["lhv"])
                                        * conversion_factor
                                    )

                                    lhv = self.fuels_specs[var]["lhv"]
                                    co2_factor = self.fuels_specs[var]["co2"]
                                    biogenic_co2_share = self.fuels_specs[var][
                                        "biogenic_share"
                                    ]

                                    fossil_co2 += (
                                        amount
                                        * lhv
                                        * co2_factor
                                        * (1 - biogenic_co2_share)
                                    )
                                    non_fossil_co2 += (
                                        amount * lhv * co2_factor * biogenic_co2_share
                                    )

                                    final_lhv += amount * lhv

                                    dataset["exchanges"].append(
                                        {
                                            "uncertainty type": 0,
                                            "amount": amount,
                                            "product": key[2],
                                            "name": key[0],
                                            "unit": key[-1],
                                            "location": key[1],
                                            "type": "technosphere",
                                        }
                                    )

                                    # log
                                    new_fuel_markets.append(
                                        [
                                            dataset["name"],
                                            dataset["location"],
                                            dataset["unit"],
                                            dataset["reference product"],
                                            fuel,
                                            key[0],
                                            key[2],
                                            key[1],
                                            key[-1],
                                            share,
                                            amount,
                                            lhv,
                                            co2_factor,
                                            biogenic_co2_share,
                                        ]
                                    )

                                    string += (
                                        f"{var.capitalize()}: {(share * 100):.1f} pct @ "
                                        f"{self.fuels_specs[var]['lhv']} MJ/kg. "
                                    )

                        string += f"Final average LHV of {final_lhv} MJ/kg."

                        if "comment" in dataset:
                            dataset["comment"] += string
                        else:
                            dataset["comment"] = string

                        # add two new fields: `fossil CO2` and `biogenic CO2`
                        dataset["fossil CO2"] = fossil_co2
                        dataset["non-fossil CO2"] = non_fossil_co2
                        dataset["LHV"] = final_lhv

                        # add fuel market to the dictionary
                        self.new_fuel_markets[
                            (dataset["name"], dataset["location"])
                        ] = {
                            "fossil CO2": fossil_co2,
                            "non-fossil CO2": non_fossil_co2,
                            "LHV": final_lhv,
                        }

                    else:
                        # World market
                        dataset["exchanges"] = [
                            e for e in dataset["exchanges"] if e["type"] == "production"
                        ]

                        final_lhv, final_fossil_co2, final_biogenic_co2 = (0, 0, 0)

                        for r in [x for x in d_act if x != "World"]:
                            total_prod_vol = (
                                self.iam_data.production_volumes.sel(
                                    region="World", variables=prod_vars
                                )
                                .interp(year=self.year)
                                .sum(dim="variables")
                                .values.item(0)
                            )

                            region_prod = (
                                self.iam_data.production_volumes.sel(
                                    region=r, variables=prod_vars
                                )
                                .interp(year=self.year)
                                .sum(dim="variables")
                                .values.item(0)
                            )

                            share = region_prod / total_prod_vol

                            dataset["exchanges"].append(
                                {
                                    "uncertainty type": 0,
                                    "amount": share,
                                    "type": "technosphere",
                                    "product": dataset["reference product"],
                                    "name": dataset["name"],
                                    "unit": dataset["unit"],
                                    "location": r,
                                }
                            )

                            lhv = self.new_fuel_markets[(dataset["name"], r)]["LHV"]
                            co2_factor = self.new_fuel_markets[(dataset["name"], r)][
                                "fossil CO2"
                            ]
                            biogenic_co2_factor = self.new_fuel_markets[
                                (dataset["name"], r)
                            ]["non-fossil CO2"]

                            final_lhv += share * lhv
                            final_fossil_co2 += share * co2_factor
                            final_biogenic_co2 += share * biogenic_co2_factor

                        # add fuel market to the dictionary
                        self.new_fuel_markets[
                            (dataset["name"], dataset["location"])
                        ] = {
                            "fossil CO2": final_fossil_co2,
                            "non-fossil CO2": final_biogenic_co2,
                            "LHV": final_lhv,
                        }

                self.database.extend(d_act.values())

        # write log to CSV
        Path(DATA_DIR / "logs").mkdir(parents=True, exist_ok=True)
        with open(
            DATA_DIR
            / f"logs/log created fuel markets {self.model} {self.scenario} {self.year}-{date.today()}.csv",
            "a",
            encoding="utf-8",
        ) as csv_file:
            writer = csv.writer(csv_file, delimiter=";", lineterminator="\n")
            for line in new_fuel_markets:
                writer.writerow(line)

        self.relink_activities_to_new_markets()

        print(f"Log of deleted fuel markets saved in {DATA_DIR / 'logs'}")
        print(f"Log of created fuel markets saved in {DATA_DIR / 'logs'}")

        print("Done!")
