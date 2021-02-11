from numpy.testing import assert_allclose


def between(bound1, bound2, value, atol=1e-5, rtol=1e-7):
    lower = min(bound1, bound2)
    upper = max(bound1, bound2)
    return lower - (rtol * lower + atol) <= value <= upper + rtol * upper + atol


def material_between(bound1, bound2, material, atol=1.0e-5, rtol=1e-7):
    return all(
        [
            between(bound1.density, bound2.density, material.density),
            between(
                bound1.primary_velocity,
                bound2.primary_velocity,
                material.primary_velocity,
            ),
            between(
                bound1.secondary_velocity,
                bound2.secondary_velocity,
                material.secondary_velocity,
            ),
            between(bound1.bulk_modulus, bound2.bulk_modulus, material.bulk_modulus),
            between(bound1.shear_modulus, bound2.shear_modulus, material.shear_modulus),
        ]
    )


def assert_similar_material(m1, m2, use_moduli=True, atol=0.0, rtol=1e-7):
    assert_allclose(m1.density, m2.density, atol=atol, rtol=rtol)
    if use_moduli:
        assert_allclose(m1.bulk_modulus, m2.bulk_modulus, atol=atol, rtol=rtol)
        assert_allclose(m1.shear_modulus, m2.shear_modulus, atol=atol, rtol=rtol)
    else:
        assert_allclose(m1.primary_velocity, m2.primary_velocity, atol=atol, rtol=rtol)
        assert_allclose(
            m1.secondary_velocity, m2.secondary_velocity, atol=atol, rtol=rtol
        )
