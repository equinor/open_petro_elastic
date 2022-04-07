from typing import Optional

import numpy as np
from pydantic import root_validator
from pydantic.dataclasses import dataclass

from .pydantic_config import PetroElasticConfig
from .vector_type import Array


@dataclass(config=PetroElasticConfig)
class Pressure:
    overburden: Array[float] = 100e6
    reference: Array[float] = 30e6
    fluid: Optional[Array[float]] = 70e6
    max_effective: Optional[Array[float]] = None
    rock: Optional[Array[float]] = None

    def __post_init_post_parse__(self):
        if self.rock is None:
            self.rock = self.fluid
        if self.fluid is None:
            self.fluid = self.rock

    @root_validator
    def at_least_one_pore_pressure(cls, values):
        if "fluid" not in values and "rock" not in values:
            raise ValueError("Need to give at least one of rock or fluid pressure.")
        return values

    @root_validator
    def positive_effective_pressure(cls, values):
        ovb = values.get("overburden")
        ref = values.get("reference")
        fld = values.get("fluid")
        rck = values.get("rock")
        if "fluid" in values:
            if np.any(ovb - fld < 0):
                raise ValueError(
                    "Overburden pressure needs to be greater than fluid pressure."
                )
        if rck is not None:
            if np.any(ovb - rck < 0):
                raise ValueError(
                    "Overburden pressure needs to be greater than rock pressure."
                )
        if np.any(ovb - ref < 0):
            raise ValueError(
                "Overburden pressure needs to be greater than reference pressure."
            )
        return values

    @property
    def effective_fluid(self):
        return self.overburden - self.fluid

    @property
    def truncated_fluid(self):
        if self.max_effective is not None:
            return np.clip(self.effective_fluid, None, self.max_effective)
        return self.effective_fluid

    @property
    def effective_reference(self):
        return self.overburden - self.reference

    @property
    def truncated_reference(self):
        if self.max_effective is not None:
            return np.clip(self.effective_reference, None, self.max_effective)
        return self.effective_reference

    @property
    def effective_rock(self):
        return self.overburden - self.rock

    @property
    def truncated_rock(self):
        if self.max_effective is not None:
            return np.clip(self.effective_rock, None, self.max_effective)
        return self.effective_rock
