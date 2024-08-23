from generators import materials, positives, ratios
from hypothesis import assume, given
from numpy.testing import assert_allclose
from predicates import between

from open_petro_elastic.material.sandstone import friable_sand, hertz_mindlin


@given(
    ratios(min_value=0.1, max_value=0.5),
    ratios(min_value=0.1, max_value=0.5),
    positives(),
    materials(),
)
def test_friable_sand_density(porosity, critical_porosity, pressure, mineral):
    assume(critical_porosity > porosity)
    sandstone = friable_sand(mineral, porosity, critical_porosity, pressure)
    assert_allclose(sandstone.density, mineral.density * (1 - porosity))


@given(
    ratios(min_value=0.1, max_value=0.5),
    ratios(min_value=0.1, max_value=0.5),
    positives(),
    materials(),
)
def test_friable_sand(porosity, critical_porosity, pressure, mineral):
    assume(critical_porosity > porosity)
    sandstone = friable_sand(mineral, porosity, critical_porosity, pressure)
    critical_mineral = hertz_mindlin(mineral, critical_porosity, pressure)
    assert between(critical_mineral.density, mineral.density, sandstone.density)
    assert between(
        critical_mineral.shear_modulus, mineral.shear_modulus, sandstone.shear_modulus
    )
    assert between(
        critical_mineral.bulk_modulus, mineral.bulk_modulus, sandstone.bulk_modulus
    )
