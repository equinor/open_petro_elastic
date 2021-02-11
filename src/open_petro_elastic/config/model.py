import abc
from typing import Optional, Union

import numpy as np
from pydantic import parse_obj_as
from pydantic.dataclasses import dataclass
from typing_extensions import Literal

from open_petro_elastic.material import Material, hashin_shtrikman_walpole
from open_petro_elastic.material.sandstone import (
    friable_sand,
    murphy_coordination_number,
    patchy_cement,
)

from .coefficients import PolyfitCoefficients, VpOverVsCoefficients
from .material_type import ArbitraryMaterial
from .pydantic_config import PetroElasticConfig
from .vector_type import Array


@dataclass(config=PetroElasticConfig)
class AbstractPressureDependencyModel(abc.ABC):
    def pressure_dependency(
        self, material, pressure, reference_pressure, porosity, mineral
    ):
        ref = self.eval(reference_pressure, porosity, mineral)
        eff = self.eval(pressure, porosity, mineral)

        if np.any(np.abs(ref) <= 1e-20):
            raise ValueError(
                "Pressure correction did not give valid dry rock, fit is zero for reference"
            )

        density_factor, primary_factor, secondary_factor = np.transpose(eff / ref)

        if self.use_moduli():
            return Material(
                density=material.density * density_factor,
                bulk_modulus=material.bulk_modulus * primary_factor,
                shear_modulus=material.shear_modulus * secondary_factor,
            )
        else:
            return Material(
                density=material.density * density_factor,
                primary_velocity=material.primary_velocity * primary_factor,
                secondary_velocity=material.secondary_velocity * secondary_factor,
            )

    @abc.abstractmethod
    def use_moduli(self):
        pass

    @abc.abstractmethod
    def eval(self, pressure, porosity, mineral):
        pass


@dataclass(config=PetroElasticConfig)
class PolyfitModel(AbstractPressureDependencyModel):
    coefficients: PolyfitCoefficients
    type: Literal["polyfit"] = "polyfit"

    def use_moduli(self):
        return self.coefficients.use_moduli

    def eval(self, pressure, porosity=None, mineral=None):
        pres = np.transpose([pressure])
        coef = np.transpose(self.coefficients.as_list)
        if coef.shape[0] == 0:
            # As it turns out, polyval does not return
            # consistent shapes for empty arrays.
            # Return the exprected 0x3 array rather than
            # the calculated 0x1 array:
            return np.transpose([[], [], []])
        else:
            return np.polyval(coef, pres)


@dataclass(config=PetroElasticConfig)
class ExpfitModel(AbstractPressureDependencyModel):
    coefficients: PolyfitCoefficients
    type: Literal["expfit"] = "expfit"

    # Todo: validate shape of coefficients

    def use_moduli(self):
        return self.coefficients.use_moduli

    def eval(self, pressure, porosity=None, mineral=None):
        coeff = np.asarray(self.coefficients.as_list)
        return np.vectorize(
            lambda p: coeff[:, 0] + coeff[:, 1] * np.exp(p / coeff[:, 2]),
            signature="()->(n)",
            otypes=[float],
        )(pressure)


@dataclass(config=PetroElasticConfig)
class LogfitModel(AbstractPressureDependencyModel):
    coefficients: PolyfitCoefficients
    type: Literal["logfit"] = "logfit"

    # Todo: validate shape of coefficients
    def use_moduli(self):
        return self.coefficients.use_moduli

    def eval(self, pressure, porosity=None, mineral=None):
        coeff = np.asarray(self.coefficients.as_list)
        return np.vectorize(
            lambda p: coeff[:, 0] + coeff[:, 1] * np.log10(p),
            signature="()->(n)",
            otypes=[float],
        )(pressure)


@dataclass(config=PetroElasticConfig)
class PowerfitModel(AbstractPressureDependencyModel):
    coefficients: VpOverVsCoefficients
    type: Literal["powerfit"] = "powerfit"

    # Todo: validate shape of coefficients

    def use_moduli(self):
        return True

    def pressure_dependency(
        self, material, pressure, reference_pressure, porosity=None, mineral=None
    ):
        ref = self.eval(reference_pressure)
        eff = self.eval(pressure)

        density_factor, primary_factor, secondary_factor = np.transpose(eff - ref)

        vp_over_vs = (
            material.primary_velocity / material.secondary_velocity + secondary_factor
        )
        bmod = material.bulk_modulus + primary_factor
        smod_reduction = vp_over_vs ** 2 - 4.0 / 3.0
        if np.any(smod_reduction <= 1e-20):
            raise ValueError("pressure correction did not result in valid dry rock")
        return Material(
            bulk_modulus=bmod,
            shear_modulus=bmod / smod_reduction,
            density=material.density + density_factor,
        )

    def eval(self, pressure, porosity=None, mineral=None):
        coeff = np.asarray(self.coefficients.as_list)
        return np.vectorize(
            lambda p: coeff[:, 0] * np.float_power(p, coeff[:, 1]),
            signature="()->(n)",
            otypes=[float],
        )(pressure)


@dataclass(config=PetroElasticConfig)
class Polyfit2dModel:
    coefficients: Optional[PolyfitCoefficients] = None
    type: Literal["polyfit"] = "polyfit"

    def use_moduli(self):
        return self.coefficients.use_moduli

    def nonporous(self, mineral):
        """
        :return: The material representing the dense framework in the
            porous dry rock.
        """
        return mineral

    def eval_material(self, pressure, porosity, mineral):
        try:
            return self.coefficients.polyfit(porosity, mineral)
        except ValueError as e:
            raise ValueError(
                f"Polynomial fit does not result in a valid dry rock material: {e}"
            )


@dataclass(config=PetroElasticConfig)
class FriableSandModel(AbstractPressureDependencyModel):
    """
    :param coordination_number: The coordination number to use for
        friable_sand. Uses
        open_petro_elastic.murphy_coordination_number(critical_porosity) by default.
    :param critical_porosity: The critical porosity of the mineral.
    :param shear_reduction: The shear_reduction parameter to use for the
        friable_sand model.
    """

    critical_porosity: Array[float] = 0.4
    shear_reduction: Array[float] = 1.0
    coordination_number: Optional[Array[float]] = None
    type: Literal["friable_sand"] = "friable_sand"

    def __post_init_post_parse__(self):
        if self.coordination_number is None:
            self.coordination_number = murphy_coordination_number(
                self.critical_porosity
            )

    def nonporous(self, mineral):
        """
        :return: The material representing the dense framework in the
            porous dry rock.
        """
        return mineral

    def use_moduli(self):
        return True

    def eval(self, pressure, porosity, mineral):
        m = self.eval_material(pressure, porosity, mineral)
        return np.transpose([m.bulk_modulus, m.shear_modulus, m.density])

    def eval_material(self, pressure, porosity, mineral):
        return friable_sand(
            mineral,
            porosity,
            self.critical_porosity,
            pressure,
            self.coordination_number,
            self.shear_reduction,
        )


DEFAULT_CEMENT = parse_obj_as(
    ArbitraryMaterial,
    {"bulk_modulus": 36.8e9, "shear_modulus": 44.0e9, "density": 2650},
)


@dataclass(config=PetroElasticConfig)
class PatchyCementModel(AbstractPressureDependencyModel):
    """
    :param cement: Material representing cement in the patchy_cement model.
    :param cement_fraction: The fraction of cement to use in the patchy_cement model.
    :param lower_bound_pressure: The lower_bound_pressure to use in the patchy_cement model.
    :param upper_bound_cement_fraction: The upper_bound_cement_fraction to use
        in the patchy_cement model.
    :param coordination_number: The coordination number to use for
        patchy_cement and friable_sand. Uses
        open_petro_elastic.murphy_coordination_number(critical_porosity) by default.
    :param critical_porosity: The critical porosity of the mineral.
    :param shear_reduction: The shear_reduction parameter to use for the
        patchy_cement and friable_sand model.
    """

    type: Literal["patchy_cement"] = "patchy_cement"
    cement: ArbitraryMaterial = DEFAULT_CEMENT
    critical_porosity: Array[float] = 0.4
    cement_fraction: Array[float] = 0.04
    lower_bound_pressure: Array[float] = 20.0e6
    upper_bound_cement_fraction: Array[float] = 0.1
    coordination_number: Optional[Array[float]] = None
    critical_porosity: Array[float] = 0.4
    shear_reduction: Array[float] = 1.0

    def __post_init_post_parse__(self):
        if self.coordination_number is None:
            self.coordination_number = murphy_coordination_number(
                self.critical_porosity
            )

    @property
    def contact_cement_porosity(self):
        """
        The amount of cement to sand at critical porosity in patchy cement.
        """
        return self.critical_porosity - self.cement_fraction

    @property
    def upper_bound_cemented_sand_porosity(self):
        """
        The lower bound on cemented sand porosity for patchy cement.
        """
        return self.critical_porosity - self.upper_bound_cement_fraction

    def nonporous(self, mineral):
        """
        :return: The material representing the dense framework in the
            porous dry rock.
        """
        return hashin_shtrikman_walpole(
            self.cement, mineral, self.upper_bound_cement_fraction
        )

    def use_moduli(self):
        return True

    def eval(self, pressure, porosity, mineral):
        m = self.eval_material(pressure, porosity, mineral)
        return np.transpose([m.bulk_modulus, m.shear_modulus, m.density])

    def eval_material(self, pressure, porosity, mineral):
        return patchy_cement(
            mineral,
            self.cement,
            porosity,
            self.contact_cement_porosity,
            self.upper_bound_cemented_sand_porosity,
            pressure,
            self.lower_bound_pressure,
            self.critical_porosity,
            self.shear_reduction,
            self.coordination_number,
        )


PressureDependencyModel = Union[
    PatchyCementModel,
    FriableSandModel,
    PowerfitModel,
    LogfitModel,
    PolyfitModel,
    ExpfitModel,
]

SandstoneModel = Union[PatchyCementModel, FriableSandModel, Polyfit2dModel]
