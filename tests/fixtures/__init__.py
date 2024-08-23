import numpy as np
import pytest

from open_petro_elastic.config import DryRock
from open_petro_elastic.material import Material
from open_petro_elastic.material.fluid import fluid_material as fluid


@pytest.fixture
def overburden_pressure():
    return 10


@pytest.fixture
def truncate_pressure():
    return 1


@pytest.fixture
def reference_pressure():
    return 1


@pytest.fixture
def max_depth():
    return 1


@pytest.fixture
def depth():
    return 1


@pytest.fixture
def reference_depth():
    return 1


@pytest.fixture
def depth_coefficients():
    return {
        "density": [[0, 1.0], [1, 0]],
        "bulk_modulus": [[0, 1.0], [1, 0]],
        "shear_modulus": [[0, 1.0], [1, 0]],
    }


@pytest.fixture
def dryrock_properties():
    return {
        "kdry": np.ones(1),
        "mydry": np.ones(1),
        "rdry": np.ones(1),
        "vpdry": np.ones(1) * ((1 + 4 / 3) ** 0.5),
        "vsdry": np.ones(1),
    }


@pytest.fixture
def dryrock_properties_material(dryrock_properties):
    return Material(
        bulk_modulus=dryrock_properties["kdry"],
        shear_modulus=dryrock_properties["mydry"],
        density=dryrock_properties["rdry"],
        primary_velocity=dryrock_properties["vpdry"],
        secondary_velocity=dryrock_properties["vsdry"],
    )


@pytest.fixture
def fluid_properties():
    return fluid(bulk_modulus=np.ones(1), density=np.ones(1))


@pytest.fixture
def fluid_input():
    return {
        "temperature": 80.0,
        "salinity": 50000.0,
        "grgas": 0.70,
        "groil": 0.70,
        "deadoildens": 850.0,
    }


@pytest.fixture
def reservoir():
    return {
        "presfluid": np.array([300.0]),
        "gor": 1.0,
        "gorgc": 1.0,
        "has_gascond": True,
        "swat": np.array([0.3]),
        "sgas": np.array([0.3]),
        "phi": 0.36,
        "presrock": 1.0,
    }


@pytest.fixture
def mineral_mix_material():
    return Material(
        bulk_modulus=36.8,
        shear_modulus=44.0,
        density=2650.0,
    )


@pytest.fixture
def polynomial_fit_pressure_coefficients():
    return {
        "density": [
            1.0,
            0.0,
        ],
        "primary_velocity": [
            20.0,
            -40.0,
        ],
        "secondary_velocity": [
            12.0,
            -20.0,
        ],
    }


@pytest.fixture
def dryrock_reference_material():
    return Material(
        bulk_modulus=1.0,
        shear_modulus=1.0,
        density=1.0,
        primary_velocity=((1 + 4 / 3) ** 0.5),
        secondary_velocity=1.0,
    )


@pytest.fixture
def velocityfit_dryrock():
    c = {
        "primary_velocity": [[2900, -1300]],
        "secondary_velocity": [[1700, -800]],
        "density": [[0.0, 0.0], [1.0, -1.0]],
    }
    return DryRock(
        porosity=np.zeros(1),
        model={"type": "polyfit", "coefficients": c},
    )


@pytest.fixture
def modulifit_dryrock():
    c = {
        "bulk_modulus": [[20.0, -40.0]],
        "shear_modulus": [[12.0, -20.0]],
        "density": [[0.0, 0.0], [1.0, -1.0]],
    }
    return DryRock(
        porosity=np.zeros(1),
        model={"type": "polyfit", "coefficients": c},
    )
