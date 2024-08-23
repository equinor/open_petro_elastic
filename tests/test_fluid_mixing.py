import numpy as np
from generators import fluids, positives, ratios
from hypothesis import assume, given
from numpy.testing import assert_allclose
from predicates import assert_similar_material

from open_petro_elastic.material import brie_fluid_mixing, wood_fluid_mixing
from open_petro_elastic.material.fluid import fluid_material


@given(fluids(), fluids(), fluids(), ratios(), ratios())
def test_bries_equation_density(water, oil, gas, water_saturation, oil_saturation):
    assume(0.0 < water_saturation + oil_saturation <= 1.0)
    gas_saturation = 1.0 - (water_saturation + oil_saturation)

    mix = brie_fluid_mixing(
        [water, oil],
        [water_saturation, oil_saturation],
        [gas],
        [gas_saturation],
        exponent=3.0,
    )

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
    mix = brie_fluid_mixing([liquid], [1.0], [gas], [0.0])
    assert_similar_material(mix, liquid, atol=1e-6)


@given(
    fluids(bulk_modulus=positives(min_value=1.1e7)),
    fluids(bulk_modulus=positives(min_value=1.1e7)),
    fluids(bulk_modulus=positives(max_value=1.0e7)),
    fluids(bulk_modulus=positives(max_value=1.0e7)),
    ratios(min_value=0.1, max_value=0.5),
    ratios(min_value=0.1, max_value=0.5),
    ratios(min_value=0.1, max_value=0.5),
    ratios(min_value=0.1, max_value=0.5),
)
def test_bries_two_liquids_two_gases(
    fluid1, fluid2, fluid3, fluid4, liquid_sat1, liquid_sat2, gas_sat1, gas_sat2
):
    assume(liquid_sat1 + liquid_sat2 > 0)
    total = liquid_sat1 + liquid_sat2 + gas_sat1 + gas_sat2
    assume(total > 0)
    liquid_sat1 /= total
    liquid_sat2 /= total
    gas_sat1 /= total
    gas_sat2 /= total
    liquids = [fluid1, fluid2]
    gases = [fluid3, fluid4]
    mix = brie_fluid_mixing(
        liquids, [liquid_sat1, liquid_sat2], gases, [gas_sat1, gas_sat2]
    )
    assert min(g.bulk_modulus for g in gases) <= mix.bulk_modulus
    assert mix.bulk_modulus <= max(liq.bulk_modulus for liq in liquids)


def test_woods_mixing_with_partial_vectorization():
    fluid1 = fluid_material(bulk_modulus=np.ones(4) * 1e7, density=np.ones(4) * 1e4)
    fluid2 = fluid_material(bulk_modulus=np.ones(4) * 1e8, density=np.ones(4) * 1e5)
    fluid3 = fluid_material(bulk_modulus=np.ones(4) * 1e9, density=np.ones(4) * 1e6)
    s0 = 0.1
    s1 = np.array([0.5, 0.6, 0.7, 0.8])
    s2 = 1.0 - s1 - s0
    saturations = [s0, s1, s2]
    mixed = wood_fluid_mixing([fluid1, fluid2, fluid3], saturations)
    assert mixed.bulk_modulus.shape == (4,)
