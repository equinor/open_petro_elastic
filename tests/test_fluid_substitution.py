import pytest
from generators import fluids, materials
from hypothesis import assume, given
from predicates import assert_similar_material

from open_petro_elastic.material import Material, fluid_substitution
from open_petro_elastic.material.sandstone import hertz_mindlin


def test_fluid_substitution(
    snapshot, mineral_mix_material, fluid_properties, dryrock_properties_material
):
    snapshot.assert_match(
        fluid_substitution(
            dryrock_properties_material, mineral_mix_material, fluid_properties, 1.0
        )
    )


@given(materials(), fluids())
def test_fluid_substitution_preserves_positives(mineral, fluid):
    assume(fluid.poisson_ratio > 0)
    assume(mineral.poisson_ratio > 0)

    assume(mineral.bulk_modulus > fluid.bulk_modulus)

    sand = hertz_mindlin(mineral, 0.4, 1.0)

    assume(mineral.bulk_modulus > sand.bulk_modulus)

    saturated = fluid_substitution(sand, mineral, fluid, 0.4)
    assert saturated.poisson_ratio > 0


@given(materials(), fluids())
def test_fluid_substitution_zero_porosity_is_mineral_identity(mineral, fluid):
    saturated = fluid_substitution(mineral, mineral, fluid, 0.0)
    assert_similar_material(saturated, mineral)


def test_fluid_substitution_errors_mineral_less_bulky_than_sand():
    mineral = Material(
        bulk_modulus=36.8,
        shear_modulus=44.0,
        density=2650.0,
    )
    sand = hertz_mindlin(mineral, 0.4, 1.0)
    fluid = Material(bulk_modulus=15.0, shear_modulus=0.0, density=1000)
    with pytest.raises(ValueError, match="moduli of materials"):
        # Note mineral and sand swapped order of arguments
        fluid_substitution(mineral, sand, fluid, 0.4)
