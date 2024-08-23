"""
The hydro-carbon gas properties are shown in figure 2 and figure 3 of
Batzle & Wang.

We can reproduce these figures:

Figure 2
========
.. plot::

    >>> import matplotlib.pyplot as plt
    >>> import numpy as np

    >>> import open_petro_elastic.material.batzle_wang as bw

    >>> fig, ax = plt.subplots()
    >>> ax.set_ylim(0.0, 0.6)

    >>> ax.set_xlabel("TEMPERATURE (°C)")
    >>> ax.set_ylabel("GAS DENSITY (g/cm³)")
    >>> xs = np.linspace(0, 350, 1000)

    >>> def plot_gas_density(gas_gravity, pressure, *args, **kwargs):
    ...     return ax.plot(
    ...         xs,
    ...         [bw.gas(x, pressure * 1e6, gas_gravity).density / 1e3 for x in xs],
    ...         *args,
    ...         **kwargs
    ...     )

    >>> plot_gas_density(0.6, 10, "k--", label="G = 0.6")
    >>> plot_gas_density(1.2, 10, "k-", label="G = 1.2")
    >>> plot_gas_density(1.2, 25, "k-")
    >>> plot_gas_density(1.2, 50, "k-")
    >>> plot_gas_density(0.6, 25, "k--")
    >>> plot_gas_density(0.6, 50, "k--")

    >>> ax.text(15.0, 0.09, "10", rotation=-20)
    >>> ax.text(15.0, 0.21, "25", rotation=-32)
    >>> ax.text(15.0, 0.265, "50", rotation=-21)

    >>> ax.text(105.0, 0.18, "10", rotation=-31)
    >>> ax.text(100.0, 0.35, "25", rotation=-32)
    >>> ax.text(100.0, 0.415, "50MPa", rotation=-21)

    >>> ax.legend()

Figure 3
========

.. plot::

    >>> import matplotlib.pyplot as plt
    >>> import numpy as np

    >>> import open_petro_elastic.material.batzle_wang as bw

    >>> fig, ax = plt.subplots()
    >>> ax.set_ylim(0.0, 650.0)
    >>> ax.set_xlabel("TEMPERATURE (°C)")
    >>> ax.set_ylabel("GAS BULK MODULUS (MPa)")
    >>> xs = np.linspace(0.0, 350, 1000)
    >>> xs_short = np.linspace(40, 350, 1000)

    >>> def plot_gas_bulk_modulus(gas_gravity, temperatures, pressure, *args, **kwargs):
    ...     return ax.plot(
    ...         temperatures,
    ...         [
    ...             bw.gas(x, pressure * 1e6, gas_gravity).bulk_modulus / 1e6
    ...             for x in temperatures
    ...         ],
    ...         *args,
    ...         **kwargs
    ...     )

    >>> plot_gas_bulk_modulus(0.6, xs, 10, "k--", label="G = 0.6")
    >>> plot_gas_bulk_modulus(1.2, xs, 10, "k-", label="G = 1.2")
    >>> plot_gas_bulk_modulus(1.2, xs, 25, "k-")
    >>> plot_gas_bulk_modulus(1.2, xs_short, 50, "k-")
    >>> plot_gas_bulk_modulus(0.6, xs, 25, "k--")
    >>> plot_gas_bulk_modulus(0.6, xs, 50, "k--")

    >>> ax.text(100.0, 25, "10", rotation=0)
    >>> ax.text(70.0, 60, "25", rotation=-2)
    >>> ax.text(45.0, 120, "50", rotation=-18)

    >>> ax.text(8.0, 22, "10", rotation=-31)
    >>> ax.text(22.0, 300, "25", rotation=-70)
    >>> ax.text(50.0, 500, "50MPa", rotation=-70)

    >>> ax.legend()

"""

from numpy import exp

from open_petro_elastic.material.conversions import celsius_to_kelvin
from open_petro_elastic.material.fluid import fluid_material as fluid

from .ideal_gas import ideal_gas


def pseudoreduced_temperature(
    absolute_temperature,
    gas_gravity,
):
    """
    calculates pseudoreduced temperature, equation 9a from Batzle & Wang [1].

    Uses relationship from

    Thomas, L. K., Hankinson, R. W., and Phillips, K. A., 1970,
    Determination of acoustic velocities for natural gas: J. Petr.
    Tech., 22, 889-892.

    :param gas_gravity: molar mass of gas relative to air molar mas.
    :param temperature: The absolute temperature of the gas in kelvin.
    :return: Pseudoreduced temperature in kelvin.
    """
    return absolute_temperature / (94.72 + 170.75 * gas_gravity)


def pseudoreduced_pressure(pressure, gas_gravity):
    """
    calculates pseudoreduced pressure, equation 9a from Batzle & Wang [1].

    Uses relationship from

    Thomas, L. K., Hankinson, R. W., and Phillips, K. A., 1970,
    Determination of acoustic velocities for natural gas: J. Petr.
    Tech., 22, 889-892.

    :param gas_gravity: molar mass of gas relative to air molar mas.
    :param pressure: Confining pressure in MPa.
    :return: Pseudoreduced pressure in MPa.
    """
    return pressure / (4.892 - 0.4048 * gas_gravity)


def compressability_factor(absolute_temperature, pressure, gas_gravity):
    """
    calculates compressability hydro-carbon gas, equation 10b and 10c from
    Batzle & Wang [1].

    :param gas_gravity: molar mass of gas relative to air molar mas.
    :param temperature: The absolute temperature of the gas in kelvin.
    :param pressure: Confining pressure in MPa.
    :return: The density of the gas in g/cc
    """
    tpr = pseudoreduced_temperature(absolute_temperature, gas_gravity)
    ppr = pseudoreduced_pressure(pressure, gas_gravity)

    return (
        (0.03 + 0.00527 * (3.5 - tpr) ** 3) * ppr
        + 0.642 * tpr
        - 0.007 * tpr**4
        - 0.52
        + 0.109
        * (3.85 - tpr) ** 2
        / exp((0.45 + 8 * (0.56 - 1 / tpr) ** 2) * ppr**1.2 / tpr)
    )


def gas_density(absolute_temperature, pressure, gas_gravity):
    """
    The density of hydro-carbon gas, using equation 10 from Batzle & Wang [1].

    :param gas_gravity: molar mass of gas relative to air molar mas.
    :param temperature: The absolute temperature of the gas in kelvin.
    :param pressure: Confining pressure in MPa.
    :return: The density of the gas in g/cc
    """
    return ideal_gas(
        absolute_temperature, pressure * 1e-6, gas_gravity
    ).density / compressability_factor(absolute_temperature, pressure, gas_gravity)


def compressability_rate_per_pseudoreduced_pressure(
    absolute_temperature, pressure, gas_gravity
):
    """
    Derivate of compressability_factor with respect to pressure.

    :param gas_gravity: molar mass of gas relative to air molar mas.
    :param temperature: The absolute temperature of the gas in kelvin.
    :param pressure: Confining pressure in MPa.
    :return: The density of the gas in g/cc
    """
    tpr = pseudoreduced_temperature(absolute_temperature, gas_gravity)
    ppr = pseudoreduced_pressure(pressure, gas_gravity)

    return (
        0.03
        + 0.00527 * (3.5 - tpr) ** 3
        - (
            0.1308
            * (0.45 + 8 * (0.56 - tpr ** (-1)) ** 2)
            * (3.85 - tpr) ** 2
            * ppr**0.2
        )
        / (exp(((0.45 + 8 * (0.56 - tpr ** (-1)) ** 2) * ppr**1.2) / tpr) * tpr)
    )


def gas_bulk_modulus(absolute_temperature, pressure, gas_gravity):
    """
    The bulk modulus of hydro-carbon gas, using equation 11 from Batzle & Wang [1].

    :param gas_gravity: molar mass of gas relative to air molar mas.
    :param temperature: The absolute temperature of the gas in kelvin.
    :param pressure: Confining pressure in MPa.
    :return: The bulk modulus of the gas in MPa.
    """
    z = compressability_factor(absolute_temperature, pressure, gas_gravity)
    dz_dppr = compressability_rate_per_pseudoreduced_pressure(
        absolute_temperature, pressure, gas_gravity
    )

    ppr = pseudoreduced_pressure(pressure, gas_gravity)

    # Equation 11b
    gamma_0 = (
        0.85
        + 5.6 / (ppr + 2)
        + 27.1 / ((ppr + 3.5) ** 2)
        - 8.7 * exp(-0.65 * (ppr + 1))
    )

    return gamma_0 * pressure / (1 - dz_dppr * ppr / z)


def gas(temperature, pressure, gas_gravity):
    """
    :param gas_gravity: molar mass of gas relative to air molar mas.
    :param pressure: Confining pressure (Pa)
    :param temperature: Temperature (celsius).
    :return: Material representing gas as defined in Batzle & Wang [1].
    """
    return fluid(
        density=gas_density(
            celsius_to_kelvin(temperature), pressure * 1e-6, gas_gravity
        ),
        bulk_modulus=gas_bulk_modulus(
            celsius_to_kelvin(temperature), pressure * 1e-6, gas_gravity
        )
        * 1e6,
    )
