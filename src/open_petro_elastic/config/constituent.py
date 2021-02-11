from typing import Optional

import numpy as np
from pydantic import validator
from pydantic.dataclasses import dataclass

from open_petro_elastic.config.fluid_model_providers import fluid_model_providers
from open_petro_elastic.material import Material

from .material_type import MaterialType
from .pydantic_config import PetroElasticConfig
from .vector_type import Array


@dataclass(config=PetroElasticConfig)
class Constituent:
    """
    A Constituent is a part of a composite material to a
    volume or mass fraction.
    :param material: The constituent material which is part of a composite material.
    :param fraction: The ratio of constituent material to composite material.
    """

    material: MaterialType
    fraction: Optional[Array[float]] = None

    def as_material(
        self,
        temperature,
        pressure,
        fluid_model_provider=fluid_model_providers["default"],
    ):
        return self.material.as_material(
            temperature, pressure, fluid_model_provider=fluid_model_provider
        )

    @validator("fraction")
    def fraction_is_ratio(cls, v):
        if v is None:
            return v
        if np.any(v < -1e-5):
            raise ValueError("Constituent fraction must be greater than zero")
        if np.any(v > 1 + 1e-5):
            raise ValueError("Constituent fraction must be less than one")
        return v

    @validator("material", pre=True)
    def wrap_a_material(cls, v):
        if isinstance(v, Material):
            return dict(
                density=v.density,
                bulk_modulus=v.bulk_modulus,
                shear_modulus=v.shear_modulus,
            )
        return v


def fix_one_fraction(constituents):
    """
    Given a list of all constituents of a composite material,
    one of which has fraction set to None, sets the value of
    that missing fraction so that the sum of fractions is 1.0.
    """
    are_none_fraction = [c for c in constituents if c.fraction is None]
    sum_rest = sum([c.fraction for c in constituents if c.fraction is not None])
    epsilon = 1e-5
    if np.any(sum_rest > 1.0 + epsilon):
        raise ValueError(
            f"Sum of constituent fractions {[c.fraction for c in constituents]} should not exceed 1.0, but is {sum_rest}"
        )
    if len(are_none_fraction) == 0:
        return
    if len(are_none_fraction) > 1:
        raise ValueError(
            f"Number of undefined fractions {[c.fraction for c in constituents]} should not be more than 1"
        )
    are_none_fraction[0].fraction = 1.0 - sum_rest
