from dataclasses import field

import numpy as np
from pydantic import conlist, validator
from pydantic.dataclasses import dataclass
from typing_extensions import Literal

from open_petro_elastic.material import brie_fluid_mixing, wood_fluid_mixing

from .constituent import Constituent, fix_one_fraction
from .fluid_model_providers import fluid_model_providers
from .material_type import BrineMaterial, GasMaterial, OilMaterial
from .pydantic_config import PetroElasticConfig
from .vector_type import Array

DEFAULT_GAS_GRAVITY = 0.7

DEFAULT_OIL = OilMaterial(
    850.0,
    100.0,
    DEFAULT_GAS_GRAVITY,
)
DEFAULT_GAS = GasMaterial(DEFAULT_GAS_GRAVITY)
DEFAULT_BRINE = BrineMaterial(50000)


@dataclass(config=PetroElasticConfig)
class Fluids:
    """
    Fluids is a composite material of several fluid materials.
    :param constituents: The list of constituent fluids.
    :param mix_method: Either "brie" or "wood" signifying whether
        mixing should be done using brie_fluid_mixing or wood_fluid_mixing.
    :param brie_constant: The exponent used in the brie function.
    :param index_of_gas: The index of gas for the fluid mixing used for brie fluid mixing.
        Has no effect on wood fluid mixing. If a constituent of type gas is given,
        defaults to index of that constituent. Otherwise, assumes first constituent is gas.
    """

    temperature: Array[float] = 80.0
    constituents: conlist(Constituent, min_items=1) = field(
        default_factory=(
            lambda: [
                Constituent(DEFAULT_OIL, 0.3),
                Constituent(DEFAULT_BRINE, 0.3),
                Constituent(DEFAULT_GAS, 0.4),
            ]
        )
    )
    mix_method: Literal["wood", "brie"] = "wood"
    brie_constant: float = 3.0
    index_of_gas: int = 0
    fluid_model: Literal[tuple(fluid_model_providers.keys())] = "default"

    def __post_init_post_parse__(self):
        fix_one_fraction(self.constituents)
        gas_constituents = [c for c in self.constituents if c.material.type == "gas"]
        if len(gas_constituents) > 0:
            self.index_of_gas = self.constituents.index(gas_constituents[0])

    @validator("constituents")
    def fractions_should_not_exceed_one(cls, v):
        fractions = np.array([c.fraction for c in v if c.fraction is not None])
        epsilon = 1e-5
        if np.any(sum(fractions) > 1.0 + epsilon):
            raise ValueError("Sum of constituent fractions should not exceed 1.0")
        return v

    def as_mixture(self, pressure):
        """
        :returns: A material representing the mixture of constituents.
        """
        fluid_model_provider = fluid_model_providers[self.fluid_model]
        materials = [
            c.material.as_material(
                self.temperature,
                pressure.fluid,
                fluid_model_provider=fluid_model_provider,
            )
            for c in self.constituents
        ]
        fractions = [c.fraction for c in self.constituents]
        # Mixing bulk modulus
        if self.mix_method == "wood":  # Method 1 - Wood's equation
            return wood_fluid_mixing(materials, fractions)
        if self.mix_method == "brie":  # Method 2 - Brie's equation with exponent 'e'
            gas_material = materials.pop(self.index_of_gas)
            fractions.pop(self.index_of_gas)
            return brie_fluid_mixing(
                materials,
                fractions,
                gas_material,
                self.brie_constant,
            )
        else:
            raise ValueError(f"Unknown mix_method: {self.mix_method}")
