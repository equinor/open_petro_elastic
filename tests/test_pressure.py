import numpy as np
import pytest
from numpy.testing import assert_allclose
from predicates import assert_similar_material

from open_petro_elastic.config import Pressure, PressureDependency
from open_petro_elastic.material import Material


def test_pressure_dependency_no_pressure_poly_fit(
    snapshot,
    dryrock_reference_material,
    overburden_pressure,
    reference_pressure,
    truncate_pressure,
    polynomial_fit_pressure_coefficients,
):
    presrock = 1.0
    pressure = Pressure(
        overburden=overburden_pressure,
        max_effective=truncate_pressure,
        reference=reference_pressure,
        rock=presrock,
    )
    snapshot.assert_match(
        PressureDependency(
            model={
                "type": "polyfit",
                "coefficients": polynomial_fit_pressure_coefficients,
            },
        ).adjust_material(dryrock_reference_material, pressure)
    )


def test_pressure_dependency_vector_pressure(dryrock_reference_material):
    pressure = Pressure(
        overburden=[20, 20],
        max_effective=[20, 20],
        reference=[5, 6],
        rock=[7, 8],
    )
    dryrock = PressureDependency(
        model={
            "type": "polyfit",
            "coefficients": {
                "bulk_modulus": [9, 10],
                "shear_modulus": [11, 12],
                "density": [13, 14],
            },
        },
    ).adjust_material(dryrock_reference_material, pressure)

    refb = dryrock_reference_material.bulk_modulus
    rockp = pressure.truncated_rock
    refp = pressure.truncated_reference
    assert_allclose(dryrock.bulk_modulus, (10 + 9 * rockp) / (10 + 9 * refp) * refb)


def test_pressure_dependency_empty_vector(dryrock_reference_material):
    pressure = Pressure(
        overburden=[],
        max_effective=[],
        reference=[],
        rock=[],
    )
    dryrock = PressureDependency(
        model={
            "type": "polyfit",
            "coefficients": {
                "bulk_modulus": [],
                "shear_modulus": [],
                "density": [],
            },
        },
    ).adjust_material(dryrock_reference_material, pressure)

    assert len(dryrock.bulk_modulus) == 0
    assert len(dryrock.density) == 0
    assert len(dryrock.shear_modulus) == 0


def test_pressure_dependency_no_pressure_power_fit(
    snapshot,
    dryrock_reference_material,
    overburden_pressure,
    reference_pressure,
    truncate_pressure,
):
    presrock = 1.0
    pressure = Pressure(
        overburden=overburden_pressure,
        max_effective=truncate_pressure,
        reference=reference_pressure,
        rock=presrock,
    )
    d = PressureDependency(
        model={
            "type": "logfit",
            "coefficients": {
                "density": np.ones(2),
                "primary_velocity": np.ones(2),
                "secondary_velocity": np.ones(2),
            },
        },
    ).adjust_material(dryrock_reference_material, pressure)
    snapshot.assert_match(d)


def test_pressure_dependency_no_pressure_exp_fit(
    snapshot,
    dryrock_reference_material,
    overburden_pressure,
    reference_pressure,
    truncate_pressure,
):
    presrock = 1.0
    pressure = Pressure(
        overburden=overburden_pressure,
        max_effective=truncate_pressure,
        reference=reference_pressure,
        rock=presrock,
    )
    d = PressureDependency(
        model={
            "type": "expfit",
            "coefficients": {
                "density": np.ones(3),
                "primary_velocity": np.ones(3),
                "secondary_velocity": np.ones(3),
            },
        },
    ).adjust_material(dryrock_reference_material, pressure)
    snapshot.assert_match(d)


def test_pressure_case_max_pressure():
    pressure = Pressure(overburden=30, fluid=10, reference=20, max_effective=11)
    assert_similar_material(
        PressureDependency(
            model={
                "type": "logfit",
                "coefficients": {
                    "density": [0.0, 1.0],
                    "bulk_modulus": [0.0, 1.0],
                    "shear_modulus": [0.0, 1.0],
                },
            },
        ).adjust_material(
            Material(bulk_modulus=1.0, shear_modulus=1.0, density=1.0), pressure
        ),
        Material(
            bulk_modulus=np.log10(11),
            shear_modulus=np.log10(11),
            density=np.log10(11),
        ),
    )


def test_pressure_invalid_reference():
    pressure = Pressure(reference=19, overburden=20, fluid=9)
    with pytest.raises(ValueError, match="zero for reference"):
        (
            PressureDependency(
                model={
                    "type": "logfit",
                    "coefficients": {
                        "density": [0.0, 1.0],
                        "bulk_modulus": [0.0, 1.0],
                        "shear_modulus": [0.0, 1.0],
                    },
                },
            ).adjust_material(
                Material(bulk_modulus=1.0, shear_modulus=1.0, density=1.0), pressure
            ),
        )


def test_pressure_invalid_powerfit_material():
    pressure = Pressure(reference=0.0, overburden=0.0, fluid=np.sqrt(3.0 / 4.0) - 10.0)
    with pytest.raises(ValueError, match="valid dry rock"):
        (
            PressureDependency(
                model={
                    "type": "powerfit",
                    "coefficients": {
                        "density": [1.0, 1.0],
                        "bulk_modulus": [1.0, 1.0],
                        "vp_over_vs": [1.0, 1.0],
                    },
                },
            ).adjust_material(
                Material(primary_velocity=10.0, secondary_velocity=1.0, density=1.0),
                pressure,
            ),
        )
