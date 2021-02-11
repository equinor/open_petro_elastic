"""
Module for calculating properties of ideal gas.
"""
from numpy import sqrt
from scipy.constants import gas_constant

from open_petro_elastic.material.fluid import fluid

AIR_WEIGHT = 28.8  # g/mol


def molecular_weight(gas_gravity):
    """
    calculates molecluar weight of a gas from gas gravity.
    :param gas_gravity: molar mass of gas relative to air molar mas.
    :return: The volume of the gas in g/mol.
    """
    return gas_gravity * AIR_WEIGHT


def molar_volume(absolute_temperature, pressure):
    """
    calculates molar volume using the ideal gas law.
    :param temperature: The absolute temperature of the gas in kelvin.
    :param pressure: Confining pressure in MPa.
    :return: The volume of the gas in cc/mol.
    """
    return gas_constant * absolute_temperature / pressure


def ideal_gas_density(absolute_temperature, pressure, gas_gravity):
    """
    calculates molar volume using the ideal gas law.
    :param gas_gravity: molar mass of gas relative to air molar mas.
    :param temperature: The absolute temperature of the gas in kelvin.
    :param pressure: Confining pressure in MPa.
    :return: The density of the gas in g/cc
    """
    return molecular_weight(gas_gravity) / molar_volume(absolute_temperature, pressure)


def ideal_gas_primary_velocity(absolute_temperature, gas_gravity):
    """
    :param gas_gravity: molar mass of gas relative to air molar mas.
    :param temperature: The absolute temperature of the gas in kelvin.
    :return: The compressional wave velocity of the gas in m/s.
    """
    return sqrt(gas_constant * absolute_temperature / molecular_weight(gas_gravity))


def ideal_gas(absolute_temperature, pressure, gas_gravity):
    """
    :param gas_gravity: molar mass of gas relative to air molar mas.
    :param temperature: The absolute temperature of the gas in kelvin.
    :param pressure: Confining pressure in Pa.
    :return: Material representing ideal gas.
    """
    return fluid(
        density=ideal_gas_density(absolute_temperature, pressure * 1e6, gas_gravity)
        * 1000,
        primary_velocity=ideal_gas_primary_velocity(absolute_temperature, gas_gravity),
    )
