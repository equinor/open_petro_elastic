"""
The velocity of water is shown in figure 12 of Batzle and Wang.  Here the
calculated values from their equations and that of the figure differ quite a
bit for 200MPa as the values in the figure are likely from experiment.  We can
still show how the calculated values fair:

Figure 12
=========

.. plot::

    >>> import matplotlib.pyplot as plt
    >>> import numpy as np

    >>> import open_petro_elastic.material.batzle_wang as bw

    >>> fig, ax = plt.subplots()
    >>> ax.set_ylim(0.5, 2.0)
    >>> ax.set_xlim(0, 350)
    >>> ax.set_xlabel("TEMPERATURE (°C)")
    >>> ax.set_ylabel("VELOCITY (km/s)")
    >>> xs = np.linspace(0.0, 350, 1000)

    >>> ax.plot(
    ...     xs,
    ...     [bw.water(x, 1.0).primary_velocity / 1000 for x in xs],
    ...     "k-",
    ...     xs,
    ...     [bw.water(x, 50).primary_velocity / 1000 for x in xs],
    ...     "k-",
    ...     xs,
    ...     [bw.water(x, 100).primary_velocity / 1000 for x in xs],
    ...     "k-",
    ... )

"""

import warnings

import numpy as np
from numpy.polynomial.polynomial import polyval2d

from open_petro_elastic.material.fluid import fluid_material as fluid


def water_density(temperature, pressure):
    """
    Density of water,, equation 27a in Batzle & Wang [1].
    :param pressure: Pressure (MPa) of oil
    :param temperature: Temperature (celsius) of oil.
    :return: Density of water in g/cc.
    """
    coefficients = [
        [1.0, 489e-6, -333e-9],
        [-8e-5, -2e-6, -2e-09],
        [-33e-7, 16e-9, 0.0],
        [1.75e-9, -13e-12, 0.0],
    ]
    return polyval2d(temperature, pressure, coefficients)


def water_primary_velocity(temperature, pressure):
    """
    Primary wave velocity of water, table 1 and equation 28 in Batzle & Wang [1].
    :param pressure: Pressure (MPa) of oil
    :param temperature: Temperature (celsius) of oil.
    :return: primary wave velocity of water in m/s.
    """
    if np.any(pressure > 100):
        warnings.warn(
            "Calculations for water velocity is not precise for\n"
            + "pressure outside [0,100]MPa"
            + f"pressure given: {pressure}MPa",
            stacklevel=1,
        )
    coefficients = [
        [1402.85, 1.524, 3.437e-3, -1.197e-5],
        [4.871, -1.11e-2, 1.739e-4, -1.628e-6],
        [-4.783e-2, 2.747e-4, -2.135e-6, 1.237e-8],
        [1.487e-4, -6.503e-7, -1.455e-8, 1.327e-10],
        [-2.197e-7, 7.987e-10, 5.23e-11, -4.614e-13],
    ]
    return polyval2d(temperature, pressure, coefficients)


def water(temperature, pressure):
    """
    :param pressure: Pressure (MPa) of oil
    :param temperature: Temperature (celsius) of oil.
    :return: Material representing brine as defined in Batzle & Wang [1].
    """
    return fluid(
        density=water_density(temperature, pressure),
        primary_velocity=water_primary_velocity(temperature, pressure),
    )
