import hypothesis.strategies as st
import pytest
from generators import materials
from hypothesis import given
from numpy.testing import assert_allclose

from open_petro_elastic.material.material import Material, polyval_material


def test_material_throws_error_if_negative():
    with pytest.raises(ValueError):
        Material(bulk_modulus=1.0, shear_modulus=1.0, density=-1.0)
    with pytest.raises(ValueError):
        Material(bulk_modulus=-1.0, shear_modulus=1.0, density=1.0)
    with pytest.raises(ValueError):
        Material(bulk_modulus=1.0, shear_modulus=-1.0, density=1.0)


@given(materials(), st.floats(allow_nan=False, allow_infinity=False))
def test_polyval_material_identity_polynomial(material, value):
    identity = [[0.0, 0.0], [1.0, 0.0]]
    new_material = polyval_material(
        value, material, identity, identity, identity, use_moduli=True
    )
    assert_allclose(material.bulk_modulus, new_material.bulk_modulus)
    assert_allclose(material.density, new_material.density)
    assert_allclose(material.shear_modulus, new_material.shear_modulus)


@given(materials(), st.floats(allow_nan=False, allow_infinity=False))
def test_polyval_material_identity_polynomial_velocity(material, value):
    identity = [[0.0, 0.0], [1.0, 0.0]]
    new_material = polyval_material(
        value, material, identity, identity, identity, use_moduli=False
    )
    assert_allclose(material.primary_velocity, new_material.primary_velocity)
    assert_allclose(material.density, new_material.density)
    assert_allclose(material.secondary_velocity, new_material.secondary_velocity)
