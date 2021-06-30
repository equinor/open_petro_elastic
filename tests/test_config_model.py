import numpy as np
import pytest
from generators import materials, positives, ratios
from hypothesis import given
from numpy.testing import assert_allclose
from predicates import assert_similar_material

from open_petro_elastic.config.model import (
    ExpfitModel,
    FriableSandModel,
    LogfitModel,
    PatchyCementModel,
    Polyfit2dModel,
    PolyfitModel,
    PowerfitModel,
)
from open_petro_elastic.material import hashin_shtrikman_walpole
from open_petro_elastic.material.sandstone import (
    friable_sand,
    murphy_coordination_number,
    patchy_cement,
)


@pytest.fixture(scope="module")
def polyfit_model():
    return PolyfitModel(
        coefficients={"bulk_modulus": [1], "shear_modulus": [1], "density": [1]},
    )


def test_polyfit_model_use_moduli(polyfit_model):
    assert polyfit_model.use_moduli()


def test_polyfit_model_eval(polyfit_model):
    assert_allclose(polyfit_model.eval(1), [1, 1, 1])


@pytest.fixture(scope="module")
def expfit_model():
    return ExpfitModel(
        coefficients={
            "bulk_modulus": [1, 0, 1],
            "shear_modulus": [1, 0, 1],
            "density": [1, 0, 1],
        },
    )


def test_expfit_model_use_moduli(expfit_model):
    assert expfit_model.use_moduli()


def test_expfit_model_eval(expfit_model):
    assert_allclose(expfit_model.eval(1), [1, 1, 1])


@pytest.fixture(scope="module")
def logfit_model():
    return LogfitModel(
        coefficients={
            "bulk_modulus": [1, 0],
            "shear_modulus": [1, 0],
            "density": [1, 0],
        },
    )


def test_logfit_model_use_moduli(logfit_model):
    assert logfit_model.use_moduli()


def test_logfit_model_eval(logfit_model):
    assert_allclose(logfit_model.eval(1), [1, 1, 1])


@pytest.fixture(scope="module")
def powerfit_model():
    return PowerfitModel(
        coefficients={
            "vp_over_vs": [1, 0],
            "bulk_modulus": [1, 0],
            "density": [1, 0],
        },
    )


def test_powerfit_model_use_moduli(powerfit_model):
    assert powerfit_model.use_moduli()


def test_powerfit_model_eval(powerfit_model):
    assert_allclose(powerfit_model.eval(1), [1, 1, 1])


@pytest.fixture(scope="module")
def polyfit2d_model():
    return Polyfit2dModel(
        coefficients={"bulk_modulus": [[1]], "shear_modulus": [[1]], "density": [[1]]},
    )


def test_polyfit2d_model_use_moduli(polyfit2d_model):
    assert polyfit2d_model.use_moduli()


@given(materials())
def test_polyfit2d_model_nonporous(polyfit2d_model, material):
    assert_similar_material(polyfit2d_model.nonporous(material), material)


def test_polyfit2d_model_eval_material(polyfit2d_model, mineral_mix_material):
    material = polyfit2d_model.eval_material(10, 0.3, mineral_mix_material)
    assert material.shear_modulus == pytest.approx(1.0)
    assert material.bulk_modulus == pytest.approx(1.0)
    assert material.density == pytest.approx(1.0)


@pytest.fixture(scope="module")
def friablesand_model():
    return FriableSandModel()


def test_friablesand_model_default_coordination_number(friablesand_model):
    assert friablesand_model.coordination_number == murphy_coordination_number(0.4)


def test_friablesand_model_use_moduli(friablesand_model):
    assert friablesand_model.use_moduli()


@given(materials())
def test_friablesand_model_nonporous(friablesand_model, material):
    assert_similar_material(friablesand_model.nonporous(material), material)


def test_friablesand_model_eval_material(friablesand_model, mineral_mix_material):
    pressure = 10e5
    porosity = 0.3
    assert_similar_material(
        friablesand_model.eval_material(pressure, porosity, mineral_mix_material),
        friable_sand(
            mineral_mix_material,
            porosity,
            0.4,
            pressure,
            murphy_coordination_number(0.4),
            1.0,
        ),
    )


@pytest.fixture(scope="module")
def patchy_cement_model():
    return PatchyCementModel()


def test_patchy_cement_model_default_contact_cement_porosity(patchy_cement_model):
    assert patchy_cement_model.contact_cement_porosity == pytest.approx(0.36)


def test_patchy_cement_model_default_upper_bound_porosity(patchy_cement_model):
    assert patchy_cement_model.upper_bound_cemented_sand_porosity == pytest.approx(0.3)


@given(materials())
def test_patchy_cement_model_nonporous(patchy_cement_model, material):
    assert_similar_material(
        patchy_cement_model.nonporous(material),
        hashin_shtrikman_walpole(
            patchy_cement_model.cement,
            material,
            patchy_cement_model.upper_bound_cement_fraction,
        ),
    )


def test_patchy_cement_model_use_moduli(patchy_cement_model):
    assert patchy_cement_model.use_moduli()


coefficient_params = {
    "coefficients": {"density": [1.0], "bulk_modulus": [1.0], "shear_modulus": [5.0]}
}


@pytest.mark.parametrize(
    "model",
    [
        PatchyCementModel(),
        FriableSandModel(),
        PolyfitModel(**coefficient_params),
        ExpfitModel(**coefficient_params),
        LogfitModel(**coefficient_params),
    ],
)
@given(materials(), ratios(min_value=0.1, max_value=0.3), positives())
def test_eval_order(model, mineral, porosity, pressure):
    model = PatchyCementModel()

    material = model.eval_material(pressure, porosity, mineral)

    assert np.all(
        model.eval(pressure, porosity, mineral)
        == np.transpose(
            [material.density, material.bulk_modulus, material.shear_modulus]
        )
    )


def test_patchy_cement_model_eval_material(patchy_cement_model, mineral_mix_material):
    pressure = 10e5
    porosity = 0.3
    assert_similar_material(
        patchy_cement_model.eval_material(pressure, porosity, mineral_mix_material),
        patchy_cement(
            mineral_mix_material,
            patchy_cement_model.cement,
            porosity,
            patchy_cement_model.contact_cement_porosity,
            patchy_cement_model.upper_bound_cemented_sand_porosity,
            pressure,
            patchy_cement_model.lower_bound_pressure,
            0.4,
            1.0,
            murphy_coordination_number(0.4),
        ),
    )
