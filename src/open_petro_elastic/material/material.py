import numpy as np
from numpy.polynomial.polynomial import polyval2d

from .seismic_velocity import (
    calculate_bulk_modulus,
    calculate_primary_velocity,
    calculate_secondary_velocity,
    calculate_shear_modulus,
)


def polyval_material(
    value,
    material,
    density_coefficients,
    primary_coefficients,
    secondary_coefficients,
    use_moduli=False,
):
    """
    :returns: New Material with properties adjusted to a 2d polynomial, i.e.

    new_density = sum( c[i][j] * old_density^i * value^j)

    :param value: One dimension of the polynomial.
    :param material: The material to adjust.
    :param density_coefficients: The coefficients of the polynomial for density
    :param primary_coefficients: The coefficients of the polynomial for the
        primary property (bulk_modulus/primary_velocity)
    :param secondary_coefficients: The coefficients of the polynomial for the
        secondary property (shear_modulus/secondary_velocity)
    :param use_moduli: Whether properties are bulk_modulus/shear_modulus or
        primary_velocity/secondary_velocity
    """

    @vectorize_material
    def fit(d, p, s, v):
        a = polyval2d(
            d,
            v,
            density_coefficients,
        )
        b = polyval2d(
            p,
            v,
            primary_coefficients,
        )
        c = polyval2d(
            s,
            v,
            secondary_coefficients,
        )
        if use_moduli:
            return Material(density=a, bulk_modulus=b, shear_modulus=c)
        else:
            return Material(density=a, primary_velocity=b, secondary_velocity=c)

    if use_moduli:
        return fit(
            material.density, material.bulk_modulus, material.shear_modulus, value
        )
    else:
        return fit(
            material.density,
            material.primary_velocity,
            material.secondary_velocity,
            value,
        )


def vectorize_material(f):
    """
    vectorize a function that returns materials so that each
    property member (density, velocity, moduli) are vectors
    of floats.
    """

    def fun(*args, **kwargs):
        vector_of_materials = np.vectorize(f, otypes=[Material])(*args, **kwargs)
        if isinstance(vector_of_materials, np.ndarray):
            if len(vector_of_materials.shape) > 0:
                return Material(
                    density=np.array([m.density for m in vector_of_materials]),
                    primary_velocity=np.array(
                        [m.primary_velocity for m in vector_of_materials]
                    ),
                    secondary_velocity=np.array(
                        [m.secondary_velocity for m in vector_of_materials]
                    ),
                )
            else:
                return vector_of_materials[()]
        else:
            return vector_of_materials

    return fun


class Material:
    """
    A isotropic material having bulk modulus,
    shear modulus, and density. These can be derived from, or
    be used to derive velocity of primary and secondary waves
    according to the usual equation for estimating primary
    and secondary wave velocities in isotropic media. See
    seismic_velocity.py
    """

    def __init__(
        self,
        bulk_modulus=None,
        shear_modulus=None,
        density=None,
        secondary_velocity=None,
        primary_velocity=None,
    ):
        """
        :param bulk_modulus: The bulk modulus of the material in Pa.
        :param shear_modulus: The shear modulus of the material in Pa.
        :param density: The density of the material in kg/m^3
        :param primary_velocity: The velocity of secondary waves
            in the material in m/s.
        :param secondary_velocity: The velocity of secondary waves
            in the material in m/s.
        :return: Mineral with given properties, properties not given
            will be calculated, see Mineral.
        """
        self.bulk_modulus = bulk_modulus
        self.shear_modulus = shear_modulus
        self.density = density
        self.primary_velocity = primary_velocity
        self.secondary_velocity = secondary_velocity
        self.set_defaults()

        if np.any(self.density < 0):
            raise ValueError("Material should not have negative density")
        if np.any(self.bulk_modulus < 0):
            raise ValueError("Material should not have negative bulk modulus")
        if np.any(self.shear_modulus < 0):
            raise ValueError("Material should not have negative shear modulus")
        if np.any(self.primary_velocity < 0):
            raise ValueError("Material should not have negative primary velocity")
        if np.any(self.secondary_velocity < 0):
            raise ValueError("Material should not have negative secondary velocity")
        if np.any(self.primary_velocity < ((4 / 3) ** 0.5) * self.secondary_velocity):
            raise ValueError(
                "primary velocity must be at least sqrt(4/3) times secondary velocity"
            )

    def set_defaults(self):
        if self.shear_modulus is None:
            if self.secondary_velocity is None or self.density is None:
                raise ValueError("Cannot calculate bulk modulus for material")
            self.shear_modulus = calculate_shear_modulus(
                self.secondary_velocity, self.density
            )
        if self.bulk_modulus is None:
            if self.primary_velocity is None or self.density is None:
                raise ValueError("Cannot calculate bulk modulus for material")
            self.bulk_modulus = calculate_bulk_modulus(
                self.primary_velocity, self.shear_modulus, self.density
            )
        if self.primary_velocity is None:
            if self.bulk_modulus is None or self.density is None:
                raise ValueError("Cannot calculate primary velocity for material")
            self.primary_velocity = calculate_primary_velocity(
                self.bulk_modulus, self.shear_modulus, self.density
            )
        if self.secondary_velocity is None:
            if self.shear_modulus is None or self.density is None:
                raise ValueError("Cannot calculate secondary velocity for material")
            self.secondary_velocity = calculate_secondary_velocity(
                self.shear_modulus, self.density
            )

    @property
    def poisson_ratio(self):
        """
        Calculates poissons ratio from bulk modulus and shear modulus using
        usual conversion formulas, see
        https://en.wikipedia.org/wiki/Template:Elastic_moduli and
        G. Mavko, T. Mukerji, J. Dvorkin. The Rock Physics Handbook. Cambridge
        University Press 2003 (paperback). ISBN 0-521-54344-4
        """
        return (3 * self.bulk_modulus - 2 * self.shear_modulus) / (
            2 * (3 * self.bulk_modulus + self.shear_modulus)
        )

    @property
    def primary_wave_modulus(self):
        """
        Calculates primary wave modulus from bulk modulus and shear modulus using
        usual conversion formulas, see
        https://en.wikipedia.org/wiki/Template:Elastic_moduli and
        G. Mavko, T. Mukerji, J. Dvorkin. The Rock Physics Handbook. Cambridge
        University Press 2003 (paperback). ISBN 0-521-54344-4
        """
        return self.bulk_modulus + 4 / 3 * self.shear_modulus

    @property
    def acoustic_impedance(self):
        return self.primary_velocity * self.density

    @property
    def shear_impedance(self):
        return self.secondary_velocity * self.density

    def __repr__(self):
        return str(self)

    def __str__(self):
        return (
            "Material("
            + (
                ""
                if self.shear_modulus is None
                else f"shear_modulus={self.shear_modulus}, "
            )
            + (
                ""
                if self.bulk_modulus is None
                else f"bulk_modulus={self.bulk_modulus}, "
            )
            + (
                ""
                if self.primary_velocity is None
                else f"primary_velocity={self.primary_velocity}, "
            )
            + (
                ""
                if self.secondary_velocity is None
                else f"secondary_velocity={self.secondary_velocity}, "
            )
            + f"density={self.density})"
        )
