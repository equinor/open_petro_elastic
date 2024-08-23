from generators import materials, positives, ratios
from hypothesis import assume, given
from numpy.testing import assert_allclose

from open_petro_elastic.material.sandstone import hertz_mindlin


@given(ratios(min_value=0.1, max_value=0.5), materials(), positives())
def test_hertz_mindlin_preserves_non_auxetic(sand_porosity, mineral, pressure):
    assume(mineral.poisson_ratio > 1e-5)
    pressurized_sand = hertz_mindlin(mineral, sand_porosity, pressure)
    assert pressurized_sand.poisson_ratio >= 0


@given(ratios(min_value=0.1, max_value=0.5), materials(), positives())
def test_hertz_mindlin_density(sand_porosity, mineral, pressure):
    pressurized_sand = hertz_mindlin(mineral, sand_porosity, pressure)
    assert_allclose(pressurized_sand.density, mineral.density * (1 - sand_porosity))


@given(ratios(min_value=0.2, max_value=0.4), materials(), positives(min_value=1.1))
def test_hertz_mindlin_normal_conditions_increases_moduli(
    sand_porosity, mineral, pressure
):
    assume(mineral.poisson_ratio > 0)
    assume(pressure > 1.0)
    sand = hertz_mindlin(mineral, sand_porosity, 1.0)
    pressurized_sand = hertz_mindlin(mineral, sand_porosity, pressure)
    assert pressurized_sand.bulk_modulus > sand.bulk_modulus
    assert pressurized_sand.shear_modulus > sand.shear_modulus
