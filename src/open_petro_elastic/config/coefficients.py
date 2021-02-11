from typing import Optional, Union

from pydantic.dataclasses import dataclass

from open_petro_elastic.material import polyval_material

from .pydantic_config import PetroElasticConfig
from .vector_type import Array


@dataclass(config=PetroElasticConfig)
class ModuliCoefficients:
    density: Optional[Array[float]] = None
    bulk_modulus: Optional[Array[float]] = None
    shear_modulus: Optional[Array[float]] = None

    @property
    def as_list(self):
        return [self.density, self.bulk_modulus, self.shear_modulus]

    @property
    def use_moduli(self):
        return True

    def polyfit(self, value, material):
        return polyval_material(
            value,
            material,
            self.density,
            self.bulk_modulus,
            self.shear_modulus,
            use_moduli=True,
        )


@dataclass(config=PetroElasticConfig)
class VelocityCoefficients:
    density: Optional[Array[float]] = None
    secondary_velocity: Optional[Array[float]] = None
    primary_velocity: Optional[Array[float]] = None

    @property
    def as_list(self):
        return [self.density, self.primary_velocity, self.secondary_velocity]

    @property
    def use_moduli(self):
        return False

    def polyfit(self, value, material):
        return polyval_material(
            value,
            material,
            self.density,
            self.primary_velocity,
            self.secondary_velocity,
            use_moduli=False,
        )


@dataclass(config=PetroElasticConfig)
class VpOverVsCoefficients:
    density: Optional[Array[float]] = None
    vp_over_vs: Optional[Array[float]] = None
    bulk_modulus: Optional[Array[float]] = None

    @property
    def as_list(self):
        return [self.density, self.bulk_modulus, self.vp_over_vs]


PolyfitCoefficients = Union[ModuliCoefficients, VelocityCoefficients]
Coefficients = Union[VpOverVsCoefficients, ModuliCoefficients, VelocityCoefficients]
