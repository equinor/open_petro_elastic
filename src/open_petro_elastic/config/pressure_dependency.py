from pydantic.dataclasses import dataclass
from typing_extensions import Literal

from .model import PressureDependencyModel
from .pydantic_config import PetroElasticConfig


@dataclass(config=PetroElasticConfig)
class PressureDependency:
    """
    Adjusts given dry rock to depth.
    :param reference_dry_rock: Dry rock to adjust at reference pressure.
    :param pressure: The pressure to adjust the dry rock to.
    :param model: One of  defaults to polyfit.
    :param coefficents: The coefficients to adjust pressure with. The resulting
        dry rock properties are set to new_property =
        f(pressure)/f(reference_pressure) * old_property, where f is one of
        f_value_log, f_value_exp and polyval2d,  depending on model choice. for powerfit
        properties are calculated by adding f(pressure) -
        f(reference_pressure) to bulk_modulus and
        primary_velocity/secondary_velocity.
    """

    model: PressureDependencyModel
    type: Literal["pressure_dependency"] = "pressure_dependency"

    def adjust_material(
        self,
        material,
        pressure,
        porosity=None,
        mineral=None,
    ):
        effective_pressure = pressure.truncated_rock
        if self.model.type == "powerfit":
            effective_pressure = pressure.rock

        return self.model.pressure_dependency(
            material,
            effective_pressure,
            pressure.truncated_reference,
            porosity,
            mineral,
        )
