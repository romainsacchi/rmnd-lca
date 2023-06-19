from premise.clean_datasets import remove_uncertainty
from premise.export import *


def test_simapro_units():
    simapro_units = get_simapro_units()
    assert simapro_units["cubic meter-year"] == "m3y"


def test_simapro_compartments():
    simapro_compartments = get_simapro_compartments()
    assert simapro_compartments["in water"] == "in water"


def test_simapro_exchange_categories():
    simapro_exchange_categories = get_simapro_category_of_exchange()
    agr = simapro_exchange_categories["agricultural trailer production"]
    assert agr["main category"] == "transport"
    assert agr["category"] == "Road\Transformation\Infrastructure"


def test_references():
    ref = load_references()
    assert (
        ref["Ethanol from maize starch"]["source"]
        == "Cozzolino, F. Life Cycle Assessment of Biofuels in EU/CH, 2018."
    )
    assert ref["Ethanol from maize starch"]["description"].startswith(
        "Production of ethanol"
    )


def test_simapro_biosphere_dict():
    s_bio = get_simapro_biosphere_dictionnary()
    assert s_bio["Propanol"] == "1-Propanol"


dummy_db = [
    {
        "name": "fake activity",
        "reference product": "fake product",
        "location": "FR",
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
                "uncertainty type": 5,
                "loc": 1.1,
                "minimum": 0.8,
                "maximum": 1.35,
            },
        ],
    },
    {
        "name": "fake activity",
        "reference product": "fake product",
        "location": "FR",
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
]


def test_remove_uncertainty():
    db = dummy_db.copy()
    db = remove_uncertainty(db)
    for ds in db:
        for exc in ds["exchanges"]:
            if "uncertainty_type" in exc:
                assert exc["uncertainty_type"] == 0


def test_check_for_duplicates():
    db = dummy_db.copy()
    db = check_for_duplicates(db)
    assert len(db) != len(dummy_db)


def test_check_amount_format():
    db = dummy_db.copy()
    db = check_amount_format(db)
    for ds in db:
        for exc in ds["exchanges"]:
            assert isinstance(exc["amount"], float)
            assert exc["amount"] == 1.0
