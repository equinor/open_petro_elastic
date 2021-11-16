import hypothesis.strategies as st
import numpy as np
import pytest
from generators import constituents
from hypothesis import given
from predicates import assert_similar_material

from open_petro_elastic.config import Fluids, Pressure
from open_petro_elastic.material import wood_fluid_mixing


@given(constituents(), constituents(fraction=st.just(None)))
def test_fluids_mixed_with_wood(oil, water):
    materials = [oil.material, water.material]
    fractions = [oil.fraction, 1.0 - oil.fraction]
    assert_similar_material(
        Fluids(constituents=[oil, water], mix_method="wood").as_mixture(Pressure()),
        wood_fluid_mixing(materials, fractions),
    )


@given(constituents(), constituents(fraction=st.just(None)))
def test_fluids_from_config(oil, water):
    oil_config = {
        "material": {
            "density": oil.material.density,
            "bulk_modulus": oil.material.bulk_modulus,
            "shear_modulus": oil.material.shear_modulus,
        },
        "fraction": oil.fraction,
    }
    water_config = {
        "material": {
            "density": water.material.density,
            "bulk_modulus": water.material.bulk_modulus,
            "shear_modulus": water.material.shear_modulus,
        },
        "fraction": water.fraction,
    }
    assert_similar_material(
        Fluids(constituents=[oil_config, water_config]).as_mixture(Pressure()),
        Fluids(constituents=[oil, water]).as_mixture(Pressure()),
    )


def test_fluids_gas_idx():
    oil_config = {
        "material": {
            "density": 1,
            "bulk_modulus": 1,
            "shear_modulus": 30,
        },
        "fraction": 0.5,
    }
    gas_config = {
        "material": {
            "type": "gas",
            "gas_gravity": 0.5,
        }
    }
    fluids = Fluids(constituents=[oil_config, gas_config], mix_method="brie")

    assert len(fluids.index_of_brie_gases()) == 1
    assert fluids.index_of_brie_gases()[0] == 1


def test_fluids_mixed_with_brie_without_gases():
    oil_config = {
        "material": {
            "density": 1,
            "bulk_modulus": 1,
            "shear_modulus": 30,
        },
        "fraction": 0.5,
    }
    brine_config = {
        "material": {"type": "brine", "salinity": 0.0},
        "fraction": 0.5,
    }
    brie_fluid = Fluids(
        constituents=[oil_config, brine_config], mix_method="brie"
    ).as_mixture(Pressure())
    wood_fluid = Fluids(
        constituents=[oil_config, brine_config], mix_method="wood"
    ).as_mixture(Pressure())
    assert_similar_material(brie_fluid, wood_fluid)


def test_fluids_mixed_with_brie_using_multiple_materials():
    brine_config1 = {
        "material": {"type": "brine", "salinity": 0.0},
        "fraction": 0.4,
    }
    brine_config2 = {
        "material": {"type": "brine", "salinity": 1e-5},
        "fraction": 0.4,
    }
    gas_config1 = {
        "material": {
            "type": "gas",
            "gas_gravity": 0.5,
        },
        "fraction": 0.1,
    }
    gas_config2 = {
        "material": {
            "type": "carbon_dioxide",
        },
        "fraction": 0.1,
    }
    mixed = Fluids(
        constituents=[brine_config1, brine_config2, gas_config1, gas_config2],
        mix_method="brie",
    )
    assert len(mixed.index_of_brie_gases()) == 2
    assert 2 in mixed.index_of_brie_gases()
    assert 3 in mixed.index_of_brie_gases()


def test_fluids_vectorized(snapshot):
    oil_config = {
        "material": {
            "density": np.ones(2),
            "bulk_modulus": np.ones(2),
            "shear_modulus": 30 * np.ones(2),
        },
        "fraction": np.array([0.5, 0.3]),
    }
    gas_config = {
        "material": {
            "type": "gas",
            "gas_gravity": 0.5 * np.ones(2),
        }
    }
    fluids = Fluids(constituents=[oil_config, gas_config], mix_method="brie")
    snapshot.assert_match(fluids.as_mixture(Pressure()))


def test_fluids_raises_fraction_exceeds():
    oil_config = {
        "material": {
            "density": 1.0,
            "bulk_modulus": 1.0,
            "shear_modulus": 30,
        },
        "fraction": 0.6,
    }
    gas_config = {
        "material": {
            "type": "gas",
            "gas_gravity": 0.5,
        },
        "fraction": 0.6,
    }
    with pytest.raises(ValueError, match="Sum of"):
        Fluids(constituents=[oil_config, gas_config])


def test_fluids_raises_fraction_not_too_strict():
    oil_config = {
        "material": {
            "density": 1.0,
            "bulk_modulus": 1.0,
            "shear_modulus": 30,
        },
        "fraction": 0.6,
    }
    gas_config = {
        "material": {
            "type": "gas",
            "gas_gravity": 0.5,
        },
        # Even though sum is above one,
        # it should be accepted as roundingerror
        "fraction": 0.4 + 1e-7,
    }

    Fluids(constituents=[oil_config, gas_config])
