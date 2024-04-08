"""
geomap.py contains the Geomap class that allows to find equivalents between
the IAM locations and ecoinvent locations.
"""

import json
from collections import defaultdict
from pathlib import Path
from typing import Any, Dict, List, Optional

import yaml
from constructive_geometries import Geomatcher

from .filesystem_constants import VARIABLES_DIR

ECO_IAM_MAPPING_FILE = VARIABLES_DIR / "missing_geography_equivalences.yaml"
TOPOLOGIES_DIR = VARIABLES_DIR / "topologies"
CONSTANTS_FILE = VARIABLES_DIR / "constants.yaml"


class Geomap:
    """
    Map ecoinvent locations to IAM regions and vice-versa.

    :ivar model: IAM model (e.g., "remind", "image")

    """

    def __init__(self, model: str) -> None:
        self.model = model
        self.geo = Geomatcher(backwards_compatible=True)
        self.constants = self.load_constants()
        self.topology = self.fetch_topology(model)
        self.additional_mappings = self.get_additional_mapping()

        self.setup_geography()

    @staticmethod
    def load_constants() -> Dict[str, Any]:
        """
        Load constants from the constants.yaml file.
        """
        with open(CONSTANTS_FILE, "r", encoding="utf-8") as stream:
            return yaml.safe_load(stream)

    @staticmethod
    def load_json(filepath: Path) -> Dict:
        """
        Load a JSON file and return its contents as a dictionary.
        """
        with open(filepath, "r", encoding="utf-8") as stream:
            return json.load(stream)

    @classmethod
    def fetch_topology(cls, model: str) -> Optional[Dict]:
        """
        Find the JSON file containing the topologies of the provided model.
        """
        topology_path = TOPOLOGIES_DIR / f"{model.lower()}-topology.json"
        if topology_path.exists():
            return cls.load_json(topology_path)

        raise FileNotFoundError(
            f"Geographical definition file for the model '{model.upper()}' not found."
        )

    @classmethod
    def get_additional_mapping(cls) -> Dict[str, dict]:
        """
        Return a dictionary with additional ecoinvent to IAM mappings.
        """
        with open(ECO_IAM_MAPPING_FILE, "r", encoding="utf-8") as stream:
            return yaml.safe_load(stream)

    def setup_geography(self) -> None:
        """
        Set up geographical definitions and additional mappings for the geomatcher.
        """
        self.geo.add_definitions(self.topology, self.model.upper(), relative=True)
        self.geo.add_definitions(
            {"World": ["GLO", "RoW"]}, self.model.upper(), relative=True
        )

        assert (self.model.upper(), "World") in self.geo.keys(), list(self.geo.keys())

        self.rev_additional_mappings = defaultdict(list)
        for ecoinvent, iam in self.additional_mappings.items():
            for iam_region in iam.values():
                self.rev_additional_mappings[iam_region].append(ecoinvent)

        self.iam_regions = [
            x[1]
            for x in list(self.geo.keys())
            if isinstance(x, tuple) and x[0] == self.model.upper()
        ]

    def iam_to_ecoinvent_location(
        self, location: str, contained: bool = True
    ) -> List[str]:
        """
        Find the corresponding ecoinvent region given an IAM region.
        :param location: name of an IAM region
        :param contained: whether only geographies that are contained within
                          the IAM region should be returned. By default, `contained` is True.
        :return: list of names of ecoinvent regions
        """
        location_tuple = (str(self.model.upper()), location)

        # Start with additional mappings that might exist
        ecoinvent_locations = []
        # first, include the missing mappings
        if location in self.rev_additional_mappings:
            ecoinvent_locations.extend(self.rev_additional_mappings[location])

        # Check if the IAM location is valid to prevent KeyError
        if location_tuple not in self.geo:
            return []

        def get_search_func(loc):
            if contained:
                return self.geo.contained(loc)
            return self.geo.intersects(loc)

        for region in get_search_func(location_tuple):
            # Skip tuple regions from unsupported models
            if (
                isinstance(region, tuple)
                and region[0].lower() in self.constants["SUPPORTED_MODELS"]
            ):
                ecoinvent_locations.append(region[1])
            elif isinstance(region, str):
                ecoinvent_locations.append(region)

        if location_tuple != (self.model.upper(), "World"):
            ecoinvent_locations = [e for e in ecoinvent_locations if e != "GLO"]

        # remove if ``location`` is in the list
        ecoinvent_locations = [e for e in ecoinvent_locations if e != location]

        return ecoinvent_locations

    def ecoinvent_to_iam_location(self, location: str) -> str:
        """
        Return an IAM region name for an ecoinvent location given.
        :param location: ecoinvent location
        :return: IAM region name
        """
        iam_locations = self.map_ecoinvent_to_iam(location)

        # Handle the case where no IAM location was found
        if not iam_locations:
            raise ValueError(
                f"No IAM location found for ecoinvent location '{location}'."
            )

        if "World" in iam_locations and len(iam_locations) > 1:
            iam_locations.remove("World")

        if len(iam_locations) == 1:
            return iam_locations[0]

        # Handle cases with multiple possible IAM regions
        return self.resolve_multiple_iam_regions(iam_locations, location)

    def map_ecoinvent_to_iam(self, location: str) -> List[str]:
        """
        Map an ecoinvent location to the corresponding IAM location(s).
        """
        # Check against the list of IAM regions and blacklist
        if location in self.iam_regions and location not in self.constants.get(
            "BLACKLIST", []
        ):
            return [location]

        # Check additional mappings
        if location in self.additional_mappings:
            mapped_location = self.additional_mappings[location].get(self.model)
            if mapped_location and mapped_location in self.iam_regions:
                return [mapped_location]

        # Find IAM regions that are within, intersect with,
        # or are contained by the ecoinvent location
        return self.find_iam_regions(location)

    def find_iam_regions(self, location: str) -> List[str]:
        """
        Find IAM regions that are within, intersect with, or are contained by an ecoinvent location.
        """
        # iterate through self.geo.within, self.geo.intersects, self.geo.contained
        # and return the first IAM region found

        for method in (self.geo.within, self.geo.intersects, self.geo.contained):
            iam_locations = [
                region[1]
                for region in method(location)
                if isinstance(region, tuple) and region[0] == self.model.upper()
            ]
            if iam_locations:
                if len(iam_locations) > 1 and "World" in iam_locations:
                    iam_locations.remove("World")
                return iam_locations
        return []

    def resolve_multiple_iam_regions(
        self, iam_locations: List[str], location: str
    ) -> str:
        """
        Resolve cases where multiple IAM regions could correspond to a single ecoinvent location.
        """
        # Order of preference to resolve IAM regions
        preferred_order = self.constants.get("PREFERRED_IAM_ORDER", [])
        for preferred in preferred_order:
            if preferred in iam_locations:
                return preferred

        raise ValueError(
            f"Multiple IAM regions found for '{location}': {iam_locations}. "
            f"None matches the preferred order: {preferred_order}."
        )
