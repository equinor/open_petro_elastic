import pytest
from predicates import assert_similar_material

from open_petro_elastic.config.constituent import Constituent, fix_one_fraction
from open_petro_elastic.material.batzle_wang import brine, gas, oil


def test_fix_one_fraction_exceeds_raises():
    with pytest.raises(ValueError, match="exceed"):
        some_gas = Constituent(material={"gas_gravity": 1.5}, fraction=0.9)
        fix_one_fraction([some_gas, some_gas])


def test_fix_two_missing_raises():
    with pytest.raises(ValueError, match="Number of undefined"):
        some_gas = Constituent(material={"gas_gravity": 1.5})
        fix_one_fraction([some_gas, some_gas])


def test_constituent_ratio_limits():
    material = {"gas_gravity": 1.5}
    Constituent(material, 0.0)
    Constituent(material, 1.0)

    # Should not be too strict on rounding
    Constituent(material, fraction=[1.0 + 1e-7])

    with pytest.raises(ValueError, match="fraction"):
        Constituent(material, fraction=[-1.0])
    with pytest.raises(ValueError, match="fraction"):
        Constituent(material, fraction=[2.0])


def test_constituent_oil_from_config():
    temperature = 80.0
    pressure = 700e5
    props = {
        "reference_density": 850.0,
        "gas_oil_ratio": 200,
        "gas_gravity": 1.5,
    }
    assert_similar_material(
        oil(temperature=temperature, pressure=pressure, **props),
        Constituent(material=props).as_material(temperature, pressure),
    )


def test_constituent_condensate_from_config():
    temperature = 80.0
    pressure = 70e5
    props = {
        "oil_reference_density": 850.0,
        "gas_oil_ratio": 200,
        "gas_gravity": 1.5,
        "type": "condensate",
    }
    with pytest.raises(NotImplementedError):
        # The default equations do not implement
        # condensate material
        Constituent(material=props).as_material(temperature, pressure)


def test_constituent_gas_from_config():
    temperature = 80.0
    pressure = 70e5
    props = {
        "gas_gravity": 1.5,
    }
    assert_similar_material(
        gas(temperature=temperature, pressure=pressure, **props),
        Constituent(material=props).as_material(temperature, pressure),
    )


def test_constituent_brine_from_config():
    temperature = 80.0
    pressure = 70e5
    props = {"salinity": 50000}
    assert_similar_material(
        brine(temperature=temperature, pressure=pressure, **props),
        Constituent(material=props).as_material(temperature, pressure),
    )
