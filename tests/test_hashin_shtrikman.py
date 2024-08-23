import pytest
from generators import materials, ratios
from hypothesis import assume, given
from predicates import assert_similar_material, between

from open_petro_elastic.material.hashin_shtrikman import (
    hashin_shtrikman_average,
    hashin_shtrikman_bound,
    hashin_shtrikman_walpole,
)


@given(materials(), materials(), ratios())
def test_hashin_shtrikman_between(material1, material2, ratio):
    bound = hashin_shtrikman_bound(material1, material2, ratio)
    assert between(material1.density, material2.density, bound.density)
    assert between(
        material1.shear_modulus, material2.shear_modulus, bound.shear_modulus
    )
    assert between(material1.bulk_modulus, material2.bulk_modulus, bound.bulk_modulus)


@given(materials(), materials())
def test_hashin_shtrikman_close1(material1, material2):
    bound = hashin_shtrikman_bound(material1, material2, 1.0)
    assert_similar_material(bound, material1)


@given(materials(), materials())
def test_hashin_shtrikman_close2(material1, material2):
    bound = hashin_shtrikman_bound(material1, material2, 0.0)
    assert_similar_material(bound, material2, atol=1e-5)


@given(materials(), materials(), ratios())
def test_hashin_shtrikman_average(material1, material2, ratio):
    bound = hashin_shtrikman_average(material1, material2, ratio)
    assert between(material1.density, material2.density, bound.density)
    assert between(
        material1.shear_modulus, material2.shear_modulus, bound.shear_modulus
    )
    assert between(material1.bulk_modulus, material2.bulk_modulus, bound.bulk_modulus)


@given(materials(), materials())
def test_hashin_shtrikman_average_close1(material1, material2):
    bound = hashin_shtrikman_average(material1, material2, 1.0)
    assert_similar_material(bound, material1, atol=1e-5)


@given(materials(), materials())
def test_hashin_shtrikman_average_close2(material1, material2):
    bound = hashin_shtrikman_average(material1, material2, 0.0)
    assert_similar_material(bound, material2, atol=1e-5)


@given(materials(), materials(), ratios())
def test_hashin_shtrikman_walpole(material1, material2, ratio):
    bound = hashin_shtrikman_average(material1, material2, ratio)
    assert between(material1.density, material2.density, bound.density)
    assert between(
        material1.shear_modulus, material2.shear_modulus, bound.shear_modulus
    )
    assert between(material1.bulk_modulus, material2.bulk_modulus, bound.bulk_modulus)


@given(materials(), materials())
def test_hashin_shtrikman_walpole_close1(material1, material2):
    bound = hashin_shtrikman_walpole(material1, material2, 1.0)
    assert_similar_material(bound, material1)


@given(materials(), materials())
def test_hashin_shtrikman_walpole_close2(material1, material2):
    bound = hashin_shtrikman_walpole(material1, material2, 0.0)
    assert_similar_material(bound, material2, atol=1e-5)


@given(materials(), materials())
def test_hashin_shtrikman_walpole_lower_less_than_upper(material1, material2):
    lower = hashin_shtrikman_walpole(material1, material2, 0.5, "lower")
    upper = hashin_shtrikman_walpole(material1, material2, 0.5, "upper")
    assert lower.bulk_modulus <= upper.bulk_modulus
    assert lower.shear_modulus <= upper.shear_modulus
    assert lower.density == upper.density


@given(materials(), materials())
def test_hashin_shtrikman_walpole_non_ratio_errors(material1, material2):
    with pytest.raises(ValueError, match="ratio"):
        hashin_shtrikman_walpole(material1, material2, 1.1)

    with pytest.raises(ValueError, match="ratio"):
        hashin_shtrikman_walpole(material1, material2, -0.1)


@given(materials(), materials())
def test_hashin_shtrikman_lower_is_less_than_upper(material1, material2):
    assume(material1.shear_modulus != pytest.approx(material2.shear_modulus))
    assume(material1.bulk_modulus != pytest.approx(material2.bulk_modulus))
    assert (
        hashin_shtrikman_walpole(
            material1,
            material2,
            0.4,
            bound="lower",
        ).bulk_modulus
        < hashin_shtrikman_walpole(
            material1,
            material2,
            0.4,
            bound="upper",
        ).bulk_modulus
    )
