import warnings
from dataclasses import field
from typing import List, Union

import numpy as np
from pydantic import root_validator
from pydantic.dataclasses import dataclass

from .depth_trend import DepthTrend
from .model import Polyfit2dModel, SandstoneModel
from .pressure_dependency import PressureDependency
from .pydantic_config import PetroElasticConfig
from .vector_type import Array

DEFAULT_BULK_COEFFICIENTS = np.array([[2900.0, -1300.0]])
DEFAULT_SHEAR_COEFFICIENTS = np.array([[1700.0, -800.0]])
DEFAULT_PRIMARY_COEFFICIENTS = np.array([[20, -40.0]])
DEFAULT_SECONDARY_COEFFICIENTS = np.array([[12.0, -20.0]])

DEFAULT_DENSITY_COEFFICIENTS = np.array([[0.0, 0.0], [1.0, -1.0]])

DEFAULT_COEFFICIENTS = {
    "density": DEFAULT_DENSITY_COEFFICIENTS,
    "bulk_modulus": DEFAULT_BULK_COEFFICIENTS,
    "shear_modulus": DEFAULT_SHEAR_COEFFICIENTS,
}


@dataclass(config=PetroElasticConfig)
class DryRock:
    """
    Dry rock is a porous material, where the pore framework is made of some
    dense material.
    :param mineral: The mineral in the dense material.
    :param porosity: The fraction of pores to dense material.
    :param model: A SandstoneModel.
    :param coefficients: The coefficients p used for polyfit model. Applied so
        that e.g. bulk_modulus of the dry rock becomes p[i,j] * porosity**i *
        mineral.bulk_modulus ** j.  see MODULI_DEFAULT_COEFFICIENTS and
        VELOCITY_DEFAULT_COEFFICIENTS for defaults.
    """

    model: SandstoneModel = Polyfit2dModel(DEFAULT_COEFFICIENTS)
    porosity: Array[float] = 0.3

    adjustments: List[Union[DepthTrend, PressureDependency]] = field(
        default_factory=list
    )

    @root_validator
    def add_default_polyfit_coefficents(cls, values):
        if "model" in values and values["model"].type == "polyfit":
            if values["model"].coefficients is None:
                values["model"].coefficients = DEFAULT_COEFFICIENTS
        return values

    @root_validator
    def warn_if_duplicate_adjustment(cls, values):
        num_pressure = 0
        num_depth = 0
        if "adjustments" in values:
            for adj in values["adjustments"]:
                if isinstance(adj, PressureDependency):
                    num_pressure += 1
                if isinstance(adj, DepthTrend):
                    num_depth += 1
        if num_pressure > 1:
            warnings.warn("Adjustments contain more than one pressure dependency")
        if num_depth > 1:
            warnings.warn("Adjustments contain more than one depth trend")
        return values

    @root_validator
    def porosity_does_not_exceed_bound_patchy(cls, values):
        epsilon = 1e-5
        if "model" in values and values["model"].type == "patchy_cement":
            cem_up_bound = (
                values["model"].critical_porosity
                - values["model"].upper_bound_cement_fraction
            )
            if "porosity" in values and np.any(
                values["porosity"] - epsilon > cem_up_bound
            ):
                raise ValueError(
                    "Porosity exceeds critical porosity minus upper"
                    + f"bound of cement fraction: {values['porosity']} vs {cem_up_bound}"
                )
        return values

    @root_validator
    def porosity_does_not_exceed_critical(cls, values):
        epsilon = 1e-5
        if (
            "porosity" in values
            and "model" in values
            and hasattr(values["model"], "critical_porosity")
            and np.any(values["porosity"] + epsilon > values["model"].critical_porosity)
        ):
            raise ValueError(
                f"Porosity exceeds critical porosity: {values['porosity']}"
            )
        return values

    def __post_init_post_parse__(self):
        if hasattr(self.model, "coefficients"):
            if self.model.coefficients.density is None:
                self.model.coefficients.density = DEFAULT_DENSITY_COEFFICIENTS
            if self.model.coefficients.use_moduli:
                if self.model.coefficients.bulk_modulus is None:
                    self.model.coefficients.bulk_modulus = DEFAULT_BULK_COEFFICIENTS
                if self.model.coefficients.shear_modulus is None:
                    self.model.coefficients.shear_modulus = DEFAULT_SHEAR_COEFFICIENTS
            else:
                if self.model.coefficients.primary_velocity is None:
                    self.model.coefficients.primary_velocity = (
                        DEFAULT_PRIMARY_COEFFICIENTS
                    )
                if self.model.coefficients.secondary_velocity is None:
                    self.model.coefficients.secondary_velocity = (
                        DEFAULT_SECONDARY_COEFFICIENTS
                    )

    def nonporous(self, mineral):
        """
        :return: The material representing the dense framework in the
            porous dry rock.
        """
        return self.model.nonporous(mineral)

    def adjust_material(self, material, pressure, mineral):
        for a in self.adjustments:
            material = a.adjust_material(material, pressure, self.porosity, mineral)
        return material

    def material(self, mineral, pressure):
        model_material = self.material_from_model(pressure.truncated_rock, mineral)
        return self.adjust_material(model_material, pressure, mineral)

    def material_from_model(self, effective_pressure, mineral):
        """
        :return: The material representing the dry rock.
        """
        return self.model.eval_material(effective_pressure, self.porosity, mineral)
