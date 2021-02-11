import pytest
from generators import materials, positives, ratios
from hypothesis import given
from numpy.testing import assert_allclose
from predicates import assert_similar_material

from open_petro_elastic.config import DepthTrend, Pressure, PressureDependency
from open_petro_elastic.config.dry_rock import (
    DEFAULT_BULK_COEFFICIENTS,
    DEFAULT_PRIMARY_COEFFICIENTS,
    DryRock,
)
from open_petro_elastic.material.sandstone import (
    friable_sand,
    murphy_coordination_number,
    patchy_cement,
)


def test_dryrock_phi_fitted_vp_vs(
    snapshot,
    overburden_pressure,
    reference_pressure,
    truncate_pressure,
    velocityfit_dryrock,
    mineral_mix_material,
):
    pressure = Pressure(
        overburden=overburden_pressure,
        reference=reference_pressure,
        max_effective=truncate_pressure,
    )
    m = velocityfit_dryrock.material(mineral_mix_material, pressure)
    dp = velocityfit_dryrock.nonporous(mineral_mix_material)
    snapshot.assert_match((m, dp))


def test_dryrock_phi_k_my(
    snapshot,
    overburden_pressure,
    reference_pressure,
    truncate_pressure,
    mineral_mix_material,
    modulifit_dryrock,
):
    pressure = Pressure(
        overburden=overburden_pressure,
        reference=reference_pressure,
        max_effective=truncate_pressure,
    )

    m = modulifit_dryrock.material(mineral_mix_material, pressure)
    dp = modulifit_dryrock.nonporous(mineral_mix_material)
    snapshot.assert_match((m, dp))


@given(
    materials(),
    ratios(min_value=0.1, max_value=0.3),
)
def test_dry_rock_polyval_positive_poisson_ratio(mineral, porosity):
    dry_rock = DryRock(porosity=porosity, model={"type": "polyfit"})
    material = dry_rock.material(mineral, Pressure())
    assert material.poisson_ratio > 0.0


@given(
    materials(
        bulk_modulus=positives(min_value=100.0),
        shear_modulus=positives(min_value=100.0),
    ),
    materials(
        bulk_modulus=positives(max_value=100.0),
        shear_modulus=positives(max_value=100.0),
    ),
    ratios(min_value=0.1, max_value=0.3),
)
def test_dry_rock_from_config_patchy_cement(mineral, cement, porosity):
    cement_config = {
        "density": cement.density,
        "bulk_modulus": cement.bulk_modulus,
        "shear_modulus": cement.shear_modulus,
    }
    dry_rock = DryRock(
        porosity=porosity, model={"cement": cement_config, "type": "patchy_cement"}
    )
    pressure = Pressure()
    material = dry_rock.material(mineral, pressure)
    pcem = patchy_cement(
        mineral,
        cement,
        porosity,
        dry_rock.model.contact_cement_porosity,
        dry_rock.model.upper_bound_cemented_sand_porosity,
        pressure.truncated_rock,
        dry_rock.model.lower_bound_pressure,
        dry_rock.model.critical_porosity,
        dry_rock.model.shear_reduction,
        dry_rock.model.coordination_number,
    )
    assert_similar_material(material, pcem)


@given(materials(), ratios(min_value=0.1, max_value=0.36))
def test_dry_rock_from_config_friable_sand(mineral, porosity):
    dry_rock = DryRock(porosity=porosity, model={"type": "friable_sand"})
    pressure = Pressure()
    material = dry_rock.material(mineral, pressure)
    f_sand = friable_sand(
        mineral,
        porosity,
        dry_rock.model.critical_porosity,
        pressure.truncated_rock,
        dry_rock.model.coordination_number,
        dry_rock.model.shear_reduction,
    )
    assert_similar_material(material, f_sand)


def test_dry_rock_too_big_porosity_raises():
    with pytest.raises(ValueError, match="Porosity exceeds"):
        DryRock(
            porosity=0.5,
            model={
                "type": "friable_sand",
                "critical_porosity": 0.4,
            },
        )


def test_dry_rock_too_big_porosity_raises_patchy():
    with pytest.raises(ValueError, match="Porosity exceeds"):
        DryRock(
            porosity=0.3,
            model={
                "critical_porosity": 0.4,
                "upper_bound_cement_fraction": 0.11,
                "type": "patchy_cement",
            },
        )


def test_dry_rock_invalid_dry_rock_raises(mineral_mix_material):
    with pytest.raises(ValueError, match="Polynomial fit"):
        dry_rock = DryRock(
            model={"coefficients": {"primary_velocity": [[-1]]}, "type": "polyfit"}
        )
        dry_rock.material(mineral_mix_material, Pressure())


def test_dry_rock_default_coefficients_modulus():
    dry_rock = DryRock(model={"coefficients": {"bulk_modulus": None}})
    assert_allclose(dry_rock.model.coefficients.bulk_modulus, DEFAULT_BULK_COEFFICIENTS)


def test_dry_rock_default_coefficients_velocity():
    dry_rock = DryRock(model={"coefficients": {"primary_velocity": None}})
    assert_allclose(
        dry_rock.model.coefficients.primary_velocity, DEFAULT_PRIMARY_COEFFICIENTS
    )


def test_dry_rock_duplicate_pressure_adjustment_warning():
    d = PressureDependency(
        model={
            "type": "polyfit",
            "coefficients": {
                "density": [0.0, 1.0],
                "bulk_modulus": [0.0, 1.0],
                "shear_modulus": [0.0, 1.0],
            },
        },
    )
    with pytest.warns(UserWarning, match="more than one pressure"):
        DryRock(
            adjustments=[d, d],
        )


def test_dry_rock_duplicate_depth_adjustment_warning():
    d = DepthTrend(
        depth=20,
        coefficients={
            "density": [0.0, 1.0],
            "bulk_modulus": [0.0, 1.0],
            "shear_modulus": [0.0, 1.0],
        },
    )
    with pytest.warns(UserWarning, match="more than one depth"):
        DryRock(
            adjustments=[d, d],
        )


def test_default_coordination_is_murphy():
    dry_rock = DryRock(model={"type": "patchy_cement"})
    assert dry_rock.model.coordination_number == murphy_coordination_number(
        dry_rock.model.critical_porosity
    )
