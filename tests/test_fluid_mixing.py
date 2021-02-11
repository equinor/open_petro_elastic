from generators import fluids, ratios
from hypothesis import assume, given
from numpy.testing import assert_allclose
from predicates import assert_similar_material

from open_petro_elastic.material import brie_fluid_mixing, wood_fluid_mixing


@given(fluids(), fluids(), fluids(), ratios(), ratios())
def test_bries_equation_density(water, oil, gas, water_saturation, oil_saturation):
    assume(0.0 < water_saturation + oil_saturation <= 1.0)
    gas_saturation = 1.0 - (water_saturation + oil_saturation)

    mix = brie_fluid_mixing([water, oil], [water_saturation, oil_saturation], gas, 3.0)

    mix_density = (
        water.density * water_saturation
        + oil.density * oil_saturation
        + gas.density * gas_saturation
    )

    assert_allclose(mix.density, mix_density, atol=0.1)


@given(fluids(), fluids(), fluids(), ratios(), ratios())
def test_woods_equation_density(water, oil, gas, water_saturation, oil_saturation):
    assume(0.0 < water_saturation + oil_saturation <= 1.0)
    gas_saturation = 1.0 - (water_saturation + oil_saturation)

    mix = wood_fluid_mixing(
        [water, oil, gas], [water_saturation, oil_saturation, gas_saturation]
    )

    mix_density = (
        water.density * water_saturation
        + oil.density * oil_saturation
        + gas.density * gas_saturation
    )

    assert_allclose(mix.density, mix_density, atol=0.1)


@given(fluids(), fluids())
def test_woods_unit_saturation_preserves_liquid(liquid, gas):
    mix = wood_fluid_mixing([liquid, gas], [1.0, 0.0])
    assert_similar_material(mix, liquid)


@given(fluids(), fluids())
def test_bries_unit_saturation_preserves_liquid(liquid, gas):

    mix = brie_fluid_mixing([liquid], [1.0], gas)
    assert_similar_material(mix, liquid, atol=1e-6)
