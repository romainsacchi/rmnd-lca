# content of test_activity_maps.py
from pathlib import Path

import pytest
from bw2data.database import DatabaseChooser

from premise.clean_datasets import DatabaseCleaner


def get_dict():
    dummy_db = {
        ("dummy_db", "6543541"): {
            "name": "fake activity",
            "reference product": "fake product",
            "location": "IAI Area, Africa",
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
        }
    }
    dummy_bio = {
        ("dummy_bio", "123"): {
            "name": "1,4-Butanediol",
            "categories": ("air", "urban air close to ground"),
            "unit": "kilogram",
        }
    }
    return dummy_db, dummy_bio


def test_presence_db():
    with pytest.raises(NameError) as wrapped_error:
        DatabaseCleaner("bla", "brightway", Path("."), version="3.8")
    assert wrapped_error.type == NameError


def test_validity_db():
    dummy_db, dummy_bio = get_dict()
    db_bio = DatabaseChooser("dummy_bio")
    db_bio.write(dummy_bio)

    db_act = DatabaseChooser("dummy_db")
    db_act.write(dummy_db)

    dbc = DatabaseCleaner("dummy_db", "brightway", Path("."), version="3.9")
    assert dbc.database[0]["name"] == "fake activity"
