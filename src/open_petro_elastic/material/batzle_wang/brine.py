import numpy as np
from numpy import sqrt
from numpy.polynomial.polynomial import polyval3d

from open_petro_elastic.material.batzle_wang.water import (
    water_density,
    water_primary_velocity,
)
from open_petro_elastic.material.fluid import fluid


def brine_density(temperature, pressure, salinity):
    """
    density of sodium chloride solutions, equation 27 in Batzle & Wang [1].
    :param salinity: Salinity of solution as weight fraction (ppm/1000000) of
        sodium chloride.
    :param pressure: Pressure (MPa) of oil
    :param temperature: Temperature (celsius) of oil.
    :return: density of solution in g/cc.
    """
    coefficients = [
        [[0.668, 3e-4], [8e-5, -13e-6], [3e-6, 0.0]],
        [[0.44, -24e-4], [-33e-4, 47e-6], [0.0, 0.0]],
    ]
    return water_density(temperature, pressure) + salinity * polyval3d(
        salinity, temperature, pressure, coefficients
    )


def brine_primary_velocity(temperature, pressure, salinity):
    """
    Primary wave velocity of sodium chloride solutions, equation 29 in Batzle & Wang [1].

    :param salinity: Salinity of solution as weight fraction (ppm/1000000) of
        sodium chloride.
    :param pressure: Pressure (MPa) of oil
    :param temperature: Temperature (celsius) of oil.
    :return: density of solution in g/cc.
    """
    coefficients = np.zeros((3, 4, 3))
    coefficients[0, 0, 0] = 1170
    coefficients[0, 1, 0] = -9.6
    coefficients[0, 2, 0] = 0.055
    coefficients[0, 3, 0] = -8.5e-5
    coefficients[0, 0, 1] = 2.6
    coefficients[0, 1, 1] = -29e-4
    coefficients[0, 0, 2] = -0.0476
    coefficients[1, 0, 0] = 780
    coefficients[1, 0, 1] = -10
    coefficients[1, 0, 2] = 0.16
    coefficients[2, 0, 0] = -820

    return water_primary_velocity(temperature, pressure) + salinity * polyval3d(
        sqrt(salinity), temperature, pressure, coefficients
    )


def brine(temperature, pressure, salinity):
    """
    :param salinity: Salinity of solution as ppm of NaCl.
    :param pressure: Pressure (Pa) of oil
    :param temperature: Temperature (celsius) of oil.
    :return: Material representing brine as defined in Batzle & Wang [1].
    """
    return fluid(
        density=brine_density(temperature, pressure * 1e-6, salinity * 1e-6) * 1000,
        primary_velocity=brine_primary_velocity(
            temperature, pressure * 1e-6, salinity * 1e-6
        ),
    )
