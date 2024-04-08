import unittest

from premise.geomap import Geomap


class TestGeomap(unittest.TestCase):
    def setUp(self):
        # This is a setup function that runs before each test.
        self.geomap = Geomap("image")

    def test_load_constants(self):
        constants = self.geomap.load_constants()
        self.assertIsInstance(constants, dict)
        self.assertIn("SUPPORTED_MODELS", constants)

    def test_fetch_topology_valid_model(self):
        topology = self.geomap.fetch_topology("image")
        self.assertIsNotNone(topology)
        self.assertIsInstance(topology, dict)

    def test_fetch_topology_invalid_model(self):
        # Assumes 'invalid_model' does not have a corresponding topology file.
        with self.assertRaises(FileNotFoundError):
            self.geomap.fetch_topology("invalid_model")

    def test_iam_to_ecoinvent_location_valid(self):
        locations = self.geomap.iam_to_ecoinvent_location("WEU")
        self.assertIsInstance(locations, list)
        self.assertGreater(len(locations), 0)

    def test_iam_to_ecoinvent_location_invalid(self):
        assert self.geomap.iam_to_ecoinvent_location("foo") == []

    def test_ecoinvent_to_iam_location_valid(self):
        iam_location = self.geomap.ecoinvent_to_iam_location("FR")
        self.assertIsInstance(iam_location, str)

    def test_ecoinvent_to_iam_location_invalid(self):
        with self.assertRaises(KeyError):
            self.geomap.ecoinvent_to_iam_location("foo")

    def test_resolve_multiple_iam_regions(self):
        self.geomap = Geomap("remind")
        iam_locations = ["EUR", "NEU"]
        location = "CH"
        resolved = self.geomap.resolve_multiple_iam_regions(iam_locations, location)
        self.assertIn(resolved, iam_locations)

    def test_map_ecoinvent_to_iam(self):
        iam_locations = self.geomap.map_ecoinvent_to_iam("IT")
        self.assertIsInstance(iam_locations, list)
        self.assertGreater(len(iam_locations), 0)
        assert len(iam_locations) == 1
        assert iam_locations[0] == "WEU"

    def test_find_iam_regions(self):
        iam_locations = self.geomap.find_iam_regions("BR")
        self.assertIsInstance(iam_locations, list)
        assert len(iam_locations) == 1, iam_locations
        assert iam_locations[0] == "BRA", iam_locations


# This allows the test to be run from the command line via `python test_geomap.py`
if __name__ == "__main__":
    unittest.main()
