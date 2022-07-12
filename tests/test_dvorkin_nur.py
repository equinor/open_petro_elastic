import numpy as np
import pytest
from generators import materials, ratios
from hypothesis import assume, given
from numpy.testing import assert_allclose

from open_petro_elastic.material import Material
from open_petro_elastic.material.sandstone import contact_cement
from open_petro_elastic.material.sandstone.contact_cement import cement_radius_ratio


@given(
    ratios(min_value=0.1, max_value=0.49),
    ratios(min_value=0.1, max_value=0.4),
    materials(),
    materials(),
)
def test_dvorkin_nur(sand_porosity, cement_fraction, cement, sand):
    assume(sand_porosity >= cement_fraction)
    assume(cement.poisson_ratio >= 0)
    assume(sand.poisson_ratio >= 0)
    cemented_sand_porosity = sand_porosity - cement_fraction
    cc = contact_cement(cement, sand, cemented_sand_porosity, sand_porosity)
    if np.allclose(cemented_sand_porosity, sand_porosity):
        assert_allclose(cc.density, sand_porosity * sand.density)


@given(ratios(min_value=0.1, max_value=0.4), ratios(min_value=0.1, max_value=0.4))
def test_radius_ratio_bounds(sand_porosity, cemented_sand_porosity):
    assume(sand_porosity > cemented_sand_porosity)

    alpha = cement_radius_ratio(cemented_sand_porosity, sand_porosity)

    assert 0.0 <= alpha <= 1.0


@given(ratios(max_value=0.9))
def test_radius_ratio_same_porosity_means_no_cement(sand_porosity):
    assert cement_radius_ratio(sand_porosity, sand_porosity) == 0.0


@given(materials(), materials())
def test_non_porosity_raises(cement, mineral):
    with pytest.raises(
        ValueError, match="cemented sand porosity must be between 0 and 1"
    ):
        contact_cement(cement, mineral, 2)


@given(materials(), materials())
def test_non_critical_porosity_raises(cement, mineral):
    with pytest.raises(ValueError, match="^sand porosity must be between 0 and 1"):
        contact_cement(cement, mineral, 0.5, -1)


@given(materials(), materials())
def test_greater_sand_porosity_than_cemented_raises(cement, mineral):
    with pytest.raises(ValueError, match="sand porosity must be greater"):
        contact_cement(cement, mineral, 0.5, 0.4)


def test_contact_cement_case():
    cement = Material(shear_modulus=4.5e9, bulk_modulus=11e9, density=2650)
    mineral = Material(shear_modulus=44e9, bulk_modulus=36e9, density=3000)
    cc = contact_cement(cement, mineral, 0.26, 0.36)
    assert cc.bulk_modulus == pytest.approx(6883899982.187891)
    assert cc.shear_modulus == pytest.approx(9010201132.579664)
