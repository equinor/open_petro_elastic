from typing import Optional, Union

from pydantic.dataclasses import dataclass
from typing_extensions import Literal

from open_petro_elastic.material import Material

from .fluid_model_providers import fluid_model_providers
from .pydantic_config import PetroElasticConfig
from .vector_type import Array


@dataclass(config=PetroElasticConfig)
class ArbitraryMaterial(Material):
    density: Array[float]
    bulk_modulus: Optional[Array[float]] = None
    shear_modulus: Optional[Array[float]] = None
    secondary_velocity: Optional[Array[float]] = None
    primary_velocity: Optional[Array[float]] = None
    type: Literal["material"] = "material"

    def __post_init_post_parse__(self):
        self.set_defaults()

    def as_material(self, temperature=None, pressure=None, fluid_model_provider=None):
        return self


@dataclass(config=PetroElasticConfig)
class OilMaterial:
    reference_density: Array[float]
    gas_oil_ratio: Array[float]
    gas_gravity: Array[float]
    type: Literal["oil"] = "oil"

    def as_material(
        self,
        temperature,
        pressure,
        fluid_model_provider=fluid_model_providers["default"],
    ):
        return fluid_model_provider.oil(
            temperature,
            pressure,
            self.reference_density,
            self.gas_oil_ratio,
            self.gas_gravity,
        )


@dataclass(config=PetroElasticConfig)
class GasMaterial:
    gas_gravity: Array[float]
    type: Literal["gas"] = "gas"

    def as_material(
        self,
        temperature,
        pressure,
        fluid_model_provider=fluid_model_providers["default"],
    ):
        return fluid_model_provider.gas(temperature, pressure, self.gas_gravity)


@dataclass(config=PetroElasticConfig)
class BrineMaterial:
    salinity: Array[float]
    type: Literal["brine"] = "brine"

    def as_material(
        self,
        temperature,
        pressure,
        fluid_model_provider=fluid_model_providers["default"],
    ):
        return fluid_model_provider.brine(temperature, pressure, self.salinity)


@dataclass(config=PetroElasticConfig)
class CondensateMaterial:
    oil_reference_density: Array[float]
    gas_oil_ratio: Array[float]
    gas_gravity: Array[float]
    type: Literal["condensate"] = "condensate"

    def as_material(
        self,
        temperature,
        pressure,
        fluid_model_provider=fluid_model_providers["default"],
    ):
        return fluid_model_provider.condensate(
            temperature,
            pressure,
            self.oil_reference_density,
            self.gas_oil_ratio,
            self.gas_gravity,
        )


MaterialType = Union[
    ArbitraryMaterial, GasMaterial, OilMaterial, BrineMaterial, CondensateMaterial
]
