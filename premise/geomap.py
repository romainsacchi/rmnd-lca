"""
geomap.py contains the Geomap class that allows to find equivalents between
the IAM locations and ecoinvent locations.
"""

import json
from typing import Dict, List, Union

import yaml
from wurst import geomatcher

from . import DATA_DIR, VARIABLES_DIR

ECO_IAM_MAPPING = VARIABLES_DIR / "missing_geography_equivalences.yaml"


def load_constants():
    """
    Load constants from the constants.yaml file.
    :return: dict
    """
    with open(VARIABLES_DIR / "constants.yaml", "r", encoding="utf-8") as stream:
        data = yaml.safe_load(stream)
    return data


constants = load_constants()


def get_additional_mapping() -> Dict[str, str]:
    """
    Return a dictionary with additional ecoinvent to IAM mappings
    """
    with open(ECO_IAM_MAPPING, "r", encoding="utf-8") as stream:
        out = yaml.safe_load(stream)

    return out


def load_json(filepath: str) -> Dict:
    """
    Load a json file.
    """
    with open(filepath, "r", encoding="utf-8") as stream:
        data = json.load(stream)
    return data


class Geomap:
    """
    Map ecoinvent locations to IAM regions and vice-versa.

    :ivar model: IAM model (e.g., "remind", "image")

    """

    def __init__(self, model: str) -> None:
        self.model = model
        self.geo = geomatcher
        self.additional_mappings = get_additional_mapping()
        self.rev_additional_mappings = {}

        if model not in ["remind", "image"]:
            if "EXTRA_TOPOLOGY" in constants:
                if model in constants["EXTRA_TOPOLOGY"]:
                    self.geo.add_definitions(
                        load_json(constants["EXTRA_TOPOLOGY"][model]),
                        self.model.upper(),
                    )
            else:
                raise ValueError(
                    f"You must provide geographical definition "
                    f"of the regions of the model {model} "
                    f"if you are not using "
                    "REMIND or IMAGE."
                )

        for key, val in self.additional_mappings.items():
            if (
                self.model.upper(),
                val[self.model],
            ) not in self.rev_additional_mappings:
                self.rev_additional_mappings[(self.model.upper(), val[self.model])] = [
                    key
                ]
            else:
                self.rev_additional_mappings[
                    (self.model.upper(), val[self.model])
                ].append(key)

        self.iam_regions = [
            x[1]
            for x in list(self.geo.keys())
            if isinstance(x, tuple) and x[0] == self.model.upper()
        ]

    def iam_to_ecoinvent_location(
        self, location: str, contained: bool = True
    ) -> Union[List[str], str]:
        """
        Find the corresponding ecoinvent region given an IAM region.
        :param location: name of a IAM region
        :param contained: whether only geographies that are contained within
        the IAM region should be returned. By default, `contained` is False,
        meaning the function also returns geographies that intersects with IAM region.
        :return: name(s) of an ecoinvent region
        """

        location = (self.model.upper(), location)

        ecoinvent_locations = []

        # first, include the missing mappings
        if location in self.rev_additional_mappings:
            ecoinvent_locations.extend(self.rev_additional_mappings[location])

        try:
            searchfunc = self.geo.contained if contained else self.geo.intersects
            for region in searchfunc(location):
                if not isinstance(region, tuple):
                    ecoinvent_locations.append(region)
                else:
                    if region[0].lower() not in constants["SUPPORTED_MODELS"]:
                        ecoinvent_locations.append(region[1])

            # Current behaviour of `intersects` is to include "GLO" in all REMIND regions.
            if location != (self.model.upper(), "World"):
                ecoinvent_locations = [e for e in ecoinvent_locations if e != "GLO"]
            return ecoinvent_locations

        except KeyError:
            print(f"Can't find location {location} using the geomatcher.")
            return ["RoW"]

    def ecoinvent_to_iam_location(self, location: str) -> str:
        """
        Return an IAM region name for an ecoinvent location given.
        Set rules in case two IAM regions are within the ecoinvent region.
        :param location: ecoinvent location
        :return: IAM region name
        """

        # First, it can be that the location is already
        # an IAM location

        list_iam_regions = [
            k[1]
            for k in list(self.geo.keys())
            if isinstance(k, tuple) and k[0].lower() == self.model.lower()
        ]

        # list of ecoinvent locations that are named
        # the same as IAM locations
        blacklist = ["ME"]

        if location in list_iam_regions and location not in blacklist:
            return location

        # Second, it can be an ecoinvent region
        # with no native mapping
        if location in self.additional_mappings:
            if self.additional_mappings[location][self.model] in self.iam_regions:
                return self.additional_mappings[location][self.model]

            # likely a case of missing "EUR" region
            raise ValueError(f"Could not find equivalent for {location}.")

        # If not, then we look for IAM regions that contain it
        try:
            iam_location = [
                r[1]
                for r in self.geo.within(location)
                if r[0] == self.model.upper() and r[1] != "World"
            ]
        except KeyError:
            print(f"Can't find location {location} using the geomatcher.")
            return "World"

        # If not, then we look for IAM regions that intersects with it
        if len(iam_location) == 0:
            iam_location = [
                r[1]
                for r in self.geo.intersects(location)
                if r[0] == self.model.upper() and r[1] != "World"
            ]

        # If not, then we look for IAM regions that are contained in it
        if len(iam_location) == 0:
            iam_location = [
                r[1]
                for r in self.geo.contained(location)
                if r[0] == self.model.upper() and r[1] != "World"
            ]

        if len(iam_location) == 0:
            print(
                f"Cannot find the IAM location for {location} from IAM model {self.model}."
            )
            return "World"

        if len(iam_location) == 1:
            return iam_location[0]

        if all(x in iam_location for x in ["NEU", "EUR"]):
            return "NEU"
        if all(x in iam_location for x in ["OAS", "JPN"]):
            return "JPN"
        if all(x in iam_location for x in ["CAZ", "USA"]):
            return "CAZ"
        if all(x in iam_location for x in ["OAS", "IND"]):
            return "IND"
        if all(x in iam_location for x in ["OAS", "CHA"]):
            return "CHA"
        if all(x in iam_location for x in ["OAS", "MEA"]):
            return "MEA"
        if all(x in iam_location for x in ["OAS", "REF"]):
            return "REF"
        if all(x in iam_location for x in ["OAS", "EUR"]):
            return "EUR"
        if all(x in iam_location for x in ["REF", "EUR"]):
            return "REF"
        if all(x in iam_location for x in ["MEA", "SSA"]):
            return "MEA"

        # more than one region is found
        print(f"More than one region found for {location}:{iam_location}")
        return "World"
