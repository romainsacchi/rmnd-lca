"""
Integrates projections regarding use of metals in the economy from:
-
"""

import uuid
from functools import lru_cache
from itertools import groupby
from pprint import pprint
from typing import Optional, Tuple

import country_converter as coco
import pandas as pd
import yaml
from _operator import itemgetter

from .export import biosphere_flows_dictionary
from .logger import create_logger
from .transformation import (
    BaseTransformation,
    Dict,
    IAMDataCollection,
    InventorySet,
    List,
    Set,
    ws,
)
from .utils import DATA_DIR

logger = create_logger("metal")


def _update_metals(scenario, version, system_model, cache=None):
    metals = Metals(
        database=scenario["database"],
        model=scenario["model"],
        pathway=scenario["pathway"],
        iam_data=scenario["iam data"],
        year=scenario["year"],
        version=version,
        system_model=system_model,
        cache=cache,
    )

    metals.create_metal_markets()
    metals.update_metals_use_in_database()
    metals.relink_datasets()
    cache = metals.cache

    return scenario, cache


def load_metals_alternative_names():
    """
    Load dataframe with alternative names for metals
    """

    filepath = DATA_DIR / "metals" / "transport_activities_mapping.yaml"

    with open(filepath, "r", encoding="utf-8") as stream:
        out = yaml.safe_load(stream)

    # this dictionary has lists as values

    # create a reversed dictionary where
    # the keys are the alternative names
    # and the values are the metals

    rev_out = {}

    for k, v in out.items():
        for i in v:
            rev_out[i] = k

    return rev_out


def load_metals_transport():
    """
    Load dataframe with metals transport
    """

    filepath = DATA_DIR / "metals" / "transport_markets_data.csv"
    df = pd.read_csv(filepath, sep=",")

    # remove rows without values under Weighted Average Distance
    df = df.loc[~df["Weighted Average Distance"].isnull()]
    # remove rows with value 0 under Weighted Average Distance
    df = df.loc[df["Weighted Average Distance"] != 0]

    df["country"] = country_short = coco.convert(df["Origin Label"], to="ISO2")

    return df


def load_mining_shares_mapping():
    """
    Load mapping between mining shares from the different sources and ecoinvent
    """

    filepath = DATA_DIR / "metals" / "mining_shares_mapping.xlsx"
    df = pd.read_excel(filepath, sheet_name="Shares_mapping")

    # replace all instances of "Year " in columns by ""
    df.columns = df.columns.str.replace("Year ", "")

    return df


def load_activities_mapping():
    """
    Load mapping for the ecoinvent exchanges to be
    updated by the new metal intensities
    """

    filepath = DATA_DIR / "metals" / "activities_mapping.xlsx"
    df = pd.read_excel(filepath, sheet_name="activities_mapping")
    return df


# Define a function to replicate rows based on the generated activity sets
def extend_dataframe(df, mapping):
    """ "
     Extend a DataFrame by duplicating rows based on a mapping dictionary.

    Parameters:
    - df (pd.DataFrame): The original DataFrame to be extended.
    - mapping (dict): A dictionary with keys corresponding to the 'technology'
                      values in the DataFrame and values that are sets of processes.
    """

    new_rows = []

    for key, processes in mapping.items():
        # Find the rows in the DataFrame where the 'technology' matches the key
        matching_rows = df[df["technology"] == key]
        # For each process in the set associated with the key, duplicate the matching rows
        for process in processes:
            temp_rows = matching_rows.copy()
            temp_rows["ecoinvent_technology"] = process
            new_rows.extend(temp_rows.to_dict("records"))
    new_df = pd.DataFrame(new_rows)

    return new_df


def get_ecoinvent_metal_factors():
    """
    Load dataframe with ecoinvent factors for metals
    and convert to xarray
    """

    filepath = DATA_DIR / "metals" / "ecoinvent_factors.csv"
    df = pd.read_csv(filepath)

    # create column "activity" as a tuple
    # of the columns "name", "product" and "location".
    df["activity"] = list(zip(df["name"], df["product"], df["location"]))
    df = df.drop(columns=["name", "product", "location"])
    df = df.melt(id_vars=["activity"], var_name="metal", value_name="value")

    # create an xarray with dimensions activity, metal and year
    ds = df.groupby(["activity", "metal"]).sum()["value"].to_xarray()

    return ds


def load_post_allocation_correction_factors():
    """
    Load yaml file with post-allocation correction factors

    """

    filepath = DATA_DIR / "metals" / "post-allocation correction" / "corrections.yaml"
    with open(filepath, "r", encoding="utf-8") as stream:
        factors = yaml.safe_load(stream)
    return factors


def fetch_mapping(filepath: str) -> dict:
    """Returns a dictionary from a YML file"""

    with open(filepath, "r", encoding="utf-8") as stream:
        mapping = yaml.safe_load(stream)
    return mapping


def rev_metals_map(mapping: dict) -> dict:
    """Returns a reversed dictionary"""

    rev_mapping = {}
    for key, val in mapping.items():
        for v in val:
            rev_mapping[v] = key
    return rev_mapping


def load_conversion_factors():
    """
    Load dataframe with conversion factors for metals
    """

    filepath = DATA_DIR / "metals" / "conversion_factors.xlsx"
    df = pd.read_excel(filepath, sheet_name="Conversion factors")
    return df


def update_exchanges(
    activity: dict,
    new_amount: float,
    new_provider: dict,
    metal: str,
) -> dict:
    """
    Update exchanges for a given activity.

    :param activity: Activity to update
    :param new_amount: Result of the calculation
    :param new_provider: Dataset of the new provider
    :param metal: Metal name
    :return: Updated activity
    """
    # fetch old amount
    old_amount = sum(
        exc["amount"]
        for exc in activity["exchanges"]
        if (exc["name"], exc.get("product"), exc["type"])
        == (
            new_provider["name"],
            new_provider["reference product"],
            "technosphere",
        )
    )

    activity["exchanges"] = [
        e
        for e in activity["exchanges"]
        if (e["name"], e.get("product"), e["type"])
        != (
            new_provider["name"],
            new_provider["reference product"],
            "technosphere",
        )
    ]

    new_exchange = {
        "uncertainty type": 0,
        "amount": new_amount,
        "product": new_provider["reference product"],
        "name": new_provider["name"],
        "unit": new_provider["unit"],
        "location": new_provider["location"],
        "type": "technosphere",
    }
    activity["exchanges"].append(new_exchange)

    # Log changes
    activity.setdefault("log parameters", {})
    activity["log parameters"].setdefault("old amount", {}).update({metal: old_amount})
    activity["log parameters"].setdefault("new amount", {}).update({metal: new_amount})

    return activity


class Metals(BaseTransformation):
    """
    Class that modifies metal demand of different technologies
    according to the Database built
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
        )

        self.version = version

        self.metals = iam_data.metals  # 1
        # Precompute the median values for each metal and origin_var for the year 2020
        self.precomputed_medians = self.metals.sel(variable="median").interp(
            year=self.year, method="nearest", kwargs={"fill_value": "extrapolate"}
        )

        self.activities_mapping = load_activities_mapping()  # 4

        self.conversion_factors = load_conversion_factors()  # 3
        # Precompute conversion factors as a dictionary for faster lookups
        self.conversion_factors_dict = self.conversion_factors.set_index("Activity")[
            "Conversion_factor"
        ].to_dict()

        inv = InventorySet(self.database, self.version)

        self.activities_metals_map: Dict[str, Set] = (
            inv.generate_metals_activities_map()
        )

        self.rev_activities_metals_map: Dict[str, str] = rev_metals_map(
            self.activities_metals_map
        )

        self.extended_dataframe = extend_dataframe(
            self.activities_mapping, self.activities_metals_map
        )
        self.extended_dataframe["final_technology"] = self.extended_dataframe.apply(
            lambda row: (
                row["demanding_process"]
                if pd.notna(row["demanding_process"]) and row["demanding_process"] != ""
                else row["ecoinvent_technology"]
            ),
            axis=1,
        )

        self.biosphere_flow_codes = biosphere_flows_dictionary(version=self.version)
        self.metals_transport = load_metals_transport()
        self.alt_names = load_metals_alternative_names()

    def update_metals_use_in_database(self):
        """
        Update the database with metals use factors.
        """

        print("Integrating metals use factors.")
        for dataset in self.database:
            if dataset["name"] in self.rev_activities_metals_map:
                origin_var = self.rev_activities_metals_map[dataset["name"]]

                self.update_metal_use(dataset, origin_var)

    @lru_cache()
    def get_metal_market_dataset(self, metal_activity_name: str):
        if pd.notna(metal_activity_name) and isinstance(metal_activity_name, str):
            metal_markets = list(
                ws.get_many(
                    self.database,
                    ws.equals("name", metal_activity_name),
                    ws.either(
                        *[ws.equals("location", loc) for loc in ["World", "GLO", "RoW"]]
                    ),
                )
            )
            metal_markets = [ds for ds in metal_markets if self.is_in_index(ds)]
            if metal_markets:
                return metal_markets[0]
            else:
                raise ws.NoResults(
                    f"Could not find dataset for metal market {metal_activity_name}"
                )
        else:
            raise ValueError(f"Invalid metal activity name: {metal_activity_name}")

    def update_metal_use(
        self,
        dataset: dict,
        technology: str,
    ) -> None:
        """
        Update metal use based on metal intensity data.
        :param dataset: dataset to adjust metal use for
        :param technology: metal intensity variable name to look up
        :return: Does not return anything. Modified in place.
        """

        # Pre-fetch relevant data to minimize DataFrame operations
        tech_rows = self.extended_dataframe.loc[
            self.extended_dataframe["ecoinvent_technology"] == dataset["name"]
        ]

        if tech_rows.empty:
            print(f"No matching rows for {dataset['name']}.")
            return

        conversion_factor = self.conversion_factors_dict.get(
            tech_rows["ecoinvent_technology"].iloc[0], None
        )
        available_metals = (
            self.precomputed_medians.sel(origin_var=technology)
            .dropna(dim="metal", how="all")["metal"]
            .values
        )

        final_technology = tech_rows["final_technology"].iloc[0]
        metal_users = ws.get_many(self.database, ws.equals("name", final_technology))

        for metal_user in metal_users:
            for metal in available_metals:
                if metal in tech_rows["Element"].values:
                    metal_row = tech_rows[tech_rows["Element"] == metal].iloc[0]
                    unit_converter = metal_row.get("unit_convertor")
                    metal_activity_name = metal_row["Activity"]

                    # Ensure that all necessary data is present
                    if (
                        pd.notna(unit_converter)
                        and pd.notna(metal_activity_name)
                        and conversion_factor
                    ):
                        median_value = self.precomputed_medians.sel(
                            metal=metal, origin_var=technology
                        ).item()
                        amount = median_value * unit_converter * conversion_factor

                        # Use a try-except block
                        # to handle the lookup of
                        # the metal market dataset once
                        try:
                            dataset_metal = self.get_metal_market_dataset(
                                metal_activity_name
                            )
                        except ws.NoResults:
                            print(f"Could not find dataset for {metal_activity_name}.")
                            continue

                        update_exchanges(metal_user, amount, dataset_metal, metal)

                    else:
                        print(
                            f"Warning: Missing data for {metal} for {dataset['name']}:"
                        )
                        if pd.isna(unit_converter):
                            print(f"- unit converter")
                        if pd.isna(metal_activity_name):
                            print(f"- activity name")
                        if not conversion_factor:
                            print(f"- conversion factor")

            self.write_log(metal_user, "updated")

    def post_allocation_correction(self):
        """
        Correct for post-allocation in the database.
        """

        factors_list = load_post_allocation_correction_factors()

        for dataset in factors_list:
            filters = [
                ws.equals("name", dataset["name"]),
                ws.equals("reference product", dataset["reference product"]),
                ws.equals("unit", dataset["unit"]),
            ]

            if "location" in dataset:
                filters.append(ws.equals("location", dataset["location"]))

            for ds in ws.get_many(
                self.database,
                *filters,
            ):
                for flow in dataset["additional flow"]:
                    ds["exchanges"].append(
                        {
                            "name": flow["name"],
                            "amount": flow["amount"],
                            "unit": flow["unit"],
                            "type": "biosphere",
                            "categories": tuple(flow["categories"].split("::")),
                            "input": (
                                "biosphere3",
                                self.biosphere_flow_codes[
                                    flow["name"],
                                    flow["categories"].split("::")[0],
                                    flow["categories"].split("::")[1],
                                    flow["unit"],
                                ],
                            ),
                        }
                    )

                if "log parameters" not in ds:
                    ds["log parameters"] = {}

                for flow in dataset["additional flow"]:
                    ds["log parameters"]["post-allocation correction"] = flow["amount"]

                self.write_log(ds, "updated")

    def create_new_mining_activity(
        self,
        name: str,
        reference_product: str,
        new_locations: dict,
        geography_mapping=None,
        shares: dict = None,
    ) -> dict:
        """
        Create a new mining activity in a new location.
        """

        geography_mapping = {
            k: v
            for k, v in geography_mapping.items()
            if not self.is_in_index(
                {"name": name, "reference product": reference_product, "location": k}
            )
        }

        # Get the original datasets
        datasets = self.fetch_proxies(
            name=name,
            ref_prod=reference_product,
            regions=new_locations.values(),
            geo_mapping=geography_mapping,
            production_variable=shares,
            exact_product_match=True,
        )

        return datasets

    @lru_cache
    def convert_long_to_short_country_name(self, country_long):
        """
        Convert long country name to short country name.
        :param country_long: Long country name
        :return: Short country name
        """
        country_short = coco.convert(country_long, to="ISO2")

        if isinstance(country_short, list):
            if country_long == "France (French Guiana)":
                country_short = "GF"
            else:
                print(f"Multiple locations found for {country_long}. Using first one.")
                country_short = country_short[0]

        if country_short not in self.geo.geo.keys():
            print(f"New location {country_short} for {country_long} not found")
            return None

        return country_short

    def get_shares(self, df: pd.DataFrame, new_locations: dict, name, ref_prod) -> dict:
        """
        Get shares of each location in the dataframe.
        :param df: Dataframe with mining shares
        :param new_locations: List of new locations
        :return: Dictionary with shares of each location
        """
        shares = {}

        # we fetch the shares for each location in df
        # and we interpolate if necessary between the columns
        # 2020 to 2030

        for long_location, short_location in new_locations.items():
            share = df.loc[df["Country"] == long_location, "2020":"2030"]

            # we interpolate depending on if self.year is between 2020 and 2030
            # otherwise, we back or forward fill

            if self.year < 2020:
                share = share.iloc[:, 0]
            elif self.year > 2030:
                share = share.iloc[:, -1]
            else:
                share = share.iloc[:, self.year - 2020]

            share = share.values[0]
            shares[(name, ref_prod, short_location)] = share

        return shares

    def get_geo_mapping(self, df: pd.DataFrame, new_locations: dict) -> dict:
        """
        Fetch the value under "Region" for each location in new_locations.
        """

        mapping = {}

        for long_location, short_location in new_locations.items():
            mapping[short_location] = df.loc[
                df["Country"] == long_location, "Region"
            ].values[0]

        return mapping

    def create_region_specific_markets(self, df: pd.DataFrame) -> List[dict]:
        new_exchanges = []

        # iterate through unique pair of process - reference product in df:
        # for each pair, create a new mining activity

        for (name, ref_prod), group in df.groupby(["Process", "Reference product"]):
            new_locations = {
                c: self.convert_long_to_short_country_name(c)
                for c in group["Country"].unique()
            }
            # remove None values
            new_locations = {k: v for k, v in new_locations.items() if v}
            # fetch shares for each location in df
            shares = self.get_shares(group, new_locations, name, ref_prod)

            geography_mapping = self.get_geo_mapping(group, new_locations)

            # if not, we create it
            datasets = self.create_new_mining_activity(
                name,
                ref_prod,
                new_locations,
                geography_mapping,
                {k[2]: v for k, v in shares.items()},
            )

            # add new datasets to database
            self.database.extend(datasets.values())
            self.add_to_index(datasets.values())

            for dataset in datasets.values():
                self.write_log(dataset, "created")

            new_exchanges.extend(
                [
                    {
                        "name": k[0],
                        "product": k[1],
                        "location": k[2],
                        "unit": "kilogram",
                        "amount": share,
                        "type": "technosphere",
                    }
                    for k, share in shares.items()
                ]
            )

        # normalize amounts to 1
        total = sum([exc["amount"] for exc in new_exchanges])
        new_exchanges = [
            {k: v for k, v in exc.items() if k != "amount"}
            | {"amount": exc["amount"] / total}
            for exc in new_exchanges
        ]

        return new_exchanges

    def create_market(self, metal, df):
        # check if market already exists
        # if so, remove it

        dataset = {
            "name": f"market for {metal[0].lower() + metal[1:]}",
            "location": "World",
            "exchanges": [
                {
                    "name": f"market for {metal[0].lower() + metal[1:]}",
                    "product": f"{metal[0].lower() + metal[1:]}",
                    "location": "World",
                    "amount": 1,
                    "type": "production",
                    "unit": "kilogram",
                }
            ],
            "reference product": f"{metal[0].lower() + metal[1:]}",
            "unit": "kilogram",
            "production amount": 1,
            "comment": "Created by premise",
            "code": str(uuid.uuid4()),
        }

        # add mining exchanges
        dataset["exchanges"].extend(self.create_region_specific_markets(df))

        # add transport exchanges
        trspt_exc = self.add_transport_to_market(dataset, metal)
        if len(trspt_exc) > 0:
            dataset["exchanges"].extend(trspt_exc)

        # filter out None
        dataset["exchanges"] = [exc for exc in dataset["exchanges"] if exc]

        # remove old market dataset
        for old_market in ws.get_many(
            self.database,
            ws.equals("name", dataset["name"]),
            ws.equals("reference product", dataset["reference product"]),
            ws.exclude(ws.equals("location", "World")),
        ):
            self.remove_from_index(old_market)
            assert (
                self.is_in_index(old_market) is False
            ), f"Market {(old_market['name'], old_market['reference product'], old_market['location'])} still in index"

        return dataset

    def add_transport_to_market(self, dataset, metal) -> list:
        excs = []

        origin_shares = {
            e["location"]: e["amount"]
            for e in dataset["exchanges"]
            if e["type"] == "technosphere"
        }

        # multiply shares with the weighted transport distance from
        # the transport dataset
        for c, share in origin_shares.items():
            if metal in self.alt_names:
                trspt_data = self.get_weighted_average_distance(c, metal)
                for i, row in trspt_data.iterrows():
                    distance = row["Weighted Average Distance"]
                    mode = row["TransportMode Label"]
                    tkm = distance / 1000 * share  # convert to tonne-kilometers x share

                    if mode == "Air":
                        name = "transport, freight, aircraft, belly-freight, long haul"
                        reference_product = "transport, freight, aircraft, long haul"
                        loc = "GLO"
                    elif mode == "Sea":
                        name = "transport, freight, sea, container ship"
                        reference_product = "transport, freight, sea, container ship"
                        loc = "GLO"
                    elif mode == "Railway":
                        name = "market group for transport, freight train"
                        reference_product = "transport, freight train"
                        loc = "GLO"
                    else:
                        name = "market for transport, freight, lorry, unspecified"
                        reference_product = "transport, freight, lorry, unspecified"
                        loc = "RoW"

                    excs.append(
                        {
                            "name": name,
                            "product": reference_product,
                            "location": loc,
                            "amount": tkm,
                            "type": "technosphere",
                            "unit": "ton kilometer",
                        }
                    )

            else:
                print(
                    f"Metal {metal} not found in alternative names. Skipping transport."
                )

        # sum up duplicates
        excs = [
            {
                "name": name,
                "product": prod,
                "location": location,
                "unit": unit,
                "type": "technosphere",
                "amount": sum([exc["amount"] for exc in exchanges]),
            }
            for (name, prod, location, unit), exchanges in groupby(
                sorted(excs, key=itemgetter("name", "product", "location", "unit")),
                key=itemgetter("name", "product", "location", "unit"),
            )
        ]

        return excs

    def get_weighted_average_distance(self, country, metal):
        return self.metals_transport.loc[
            (self.metals_transport["country"] == country)
            & (self.metals_transport["Metal"] == self.alt_names[metal]),
            ["TransportMode Label", "Weighted Average Distance"],
        ]

    def create_metal_markets(self):
        self.post_allocation_correction()

        print("Creating metal markets")

        dataframe = load_mining_shares_mapping()
        dataframe = dataframe.loc[dataframe["Work done"] == "Yes"]
        dataframe = dataframe.loc[~dataframe["Country"].isnull()]
        dataframe_shares = dataframe

        for metal in dataframe_shares["Metal"].unique():
            print(f"... for {metal}.")
            df_metal = dataframe.loc[dataframe["Metal"] == metal]
            dataset = self.create_market(metal, df_metal)

            self.database.append(dataset)
            self.add_to_index(dataset)
            self.write_log(dataset, "created")

        # filter dataframe_parent to only keep rows
        # which have a region and no values under columns from 2020 to 2030

        dataframe_parent = dataframe.loc[
            dataframe["Region"].notnull() & dataframe["2020"].isnull()
        ]

        print("Creating additional mining processes")
        for metal in dataframe_parent["Metal"].unique():
            df_metal = dataframe_parent.loc[dataframe["Metal"] == metal]

            for (name, ref_prod), group in df_metal.groupby(
                ["Process", "Reference product"]
            ):
                print(f"...... for {name} - {ref_prod}.")
                new_locations = {
                    c: self.convert_long_to_short_country_name(c)
                    for c in group["Country"].unique()
                }
                # remove None
                new_locations = {
                    k: v for k, v in new_locations.items() if v is not None
                }

                geography_mapping = self.get_geo_mapping(group, new_locations)

                # if not, we create it
                datasets = self.create_new_mining_activity(
                    name, ref_prod, new_locations, geography_mapping=geography_mapping
                )

                self.database.extend(datasets.values())
                for dataset in datasets.values():
                    self.add_to_index(dataset)
                    self.write_log(dataset, "created")

    def write_log(self, dataset, status="created"):
        """
        Write log file.
        """

        txt = (
            f"{status}|{self.model}|{self.scenario}|{self.year}|"
            f"{dataset['name']}|{dataset['reference product']}|{dataset['location']}|"
            f"{dataset.get('log parameters', {}).get('post-allocation correction', '')}|"
            f"{dataset.get('log parameters', {}).get('old amount', '')}|"
            f"{dataset.get('log parameters', {}).get('new amount', '')}"
        )

        logger.info(txt)
