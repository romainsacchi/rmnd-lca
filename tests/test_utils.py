from unittest.mock import call, patch

from premise import __version__
from premise.utils import *


def test_ei_db_label():
    model = "remind"
    pathway = "SSP2-Base"
    year = 2012
    version = "3.9"
    system_model = "cutoff"
    assert (
        eidb_label(model, pathway, year, version, system_model)
        == f"ecoinvent_{system_model}_{version}_{model}_{pathway}_{year}"
    )


def test_crops_properties():
    crop_props = get_crops_properties()
    assert type(crop_props) == dict
    assert crop_props["sugar"]["crop_type"]["image"]["temperate"] == "sugarbeet"


def test_fuels_properties():
    fuels_props = get_fuel_properties()
    assert type(fuels_props) == dict
    assert fuels_props["bioethanol, wood, with CCS"]["lhv"] == 26.5


def test_eff_solar_PV():
    eff_PV = get_efficiency_solar_photovoltaics()
    assert type(eff_PV) == xr.DataArray
    assert eff_PV.sel(technology="multi-Si", year=2010) == 0.14
    assert "moni-Si" not in eff_PV.technology.values


def test_default_location():
    dummy_db = [
        {
            "name": "fake activity",
            "reference product": "fake product",
            "location": None,
            "unit": "kilogram",
            "exchanges": [
                {
                    "name": "fake activity",
                    "product": "fake product",
                    "amount": 1,
                    "type": "production",
                    "unit": "kilogram",
                    "input": ("dummy_db", "6543541"),
                },
                {
                    "name": "1,4-Butanediol",
                    "categories": ("air", "urban air close to ground"),
                    "amount": 1,
                    "type": "biosphere",
                    "unit": "kilogram",
                    "input": ("dummy_bio", "123"),
                },
            ],
        },
        {
            "name": "fake activity 2",
            "reference product": "fake product",
            "location": "FR",
            "unit": "kilogram",
            "exchanges": [
                {
                    "name": "fake activity 2",
                    "product": "fake product 2",
                    "amount": 1,
                    "type": "production",
                    "unit": "kilogram",
                    "input": ("dummy_db", "6543541"),
                },
                {
                    "name": "1,4-Butanediol",
                    "categories": ("air", "urban air close to ground"),
                    "amount": 1,
                    "type": "biosphere",
                    "unit": "kilogram",
                    "input": ("dummy_bio", "123"),
                },
            ],
        },
    ]

    corrected_db = default_global_location(dummy_db)

    assert not any(ds["location"] is None for ds in corrected_db)


@patch("builtins.print")
def test_print_version(mocked_print):
    print_version()
    assert mocked_print.mock_calls == [call(f"premise v.{__version__}")]
