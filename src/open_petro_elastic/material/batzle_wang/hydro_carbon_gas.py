from numpy import exp

from open_petro_elastic.material.conversions import celsius_to_kelvin
from open_petro_elastic.material.fluid import fluid

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
        - 0.007 * tpr ** 4
        - 0.52
        + 0.109
        * (3.85 - tpr) ** 2
        / exp((0.45 + 8 * (0.56 - 1 / tpr) ** 2) * ppr ** 1.2 / tpr)
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
            * ppr ** 0.2
        )
        / (exp(((0.45 + 8 * (0.56 - tpr ** (-1)) ** 2) * ppr ** 1.2) / tpr) * tpr)
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

    return gamma_0 * pressure / ((1 - dz_dppr * ppr / z))


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
