import numpy as np

try:
    from pydantic.v1 import ValidationError, parse_obj_as
except ImportError:
    from pydantic import ValidationError, parse_obj_as

import pytest

from open_petro_elastic.config.coefficients import Coefficients


def test_read_coefficients_use_moduli_true_order():
    c = parse_obj_as(
        Coefficients, {"bulk_modulus": [2], "shear_modulus": [3], "density": [1]}
    )
    assert c.use_moduli
    assert np.all(c.as_list == [[1], [2], [3]])


def test_read_coefficients_use_velocity_false():
    c = parse_obj_as(
        Coefficients,
        {"primary_velocity": [2], "secondary_velocity": [3], "density": [1]},
    )
    assert not c.use_moduli
    assert np.all(c.as_list == [[1], [2], [3]])


def test_read_coefficients_use_velocity_mix_throws():
    with pytest.raises(ValidationError):
        parse_obj_as(
            Coefficients, {"some other string": [], "shear_modulus": [], "density": []}
        )
