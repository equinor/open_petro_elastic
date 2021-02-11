import warnings

import numpy as np
from open_petro_elastic.material.fluid import fluid


def oil_bubble_point(density, gas_oil_ratio, gas_gravity, temperature):
    """
    Reservoir oils include some natural gas in solution. The pressure at which
    this natural gas begins to come out of solution and form bubbles is known
    as the bubble point pressure. See https://petrowiki.org/Oil_bubblepoint_pressure
    Based on the correlation from:
    Standing, M. B. "A pressure-volume-temperature correlation for mixtures of
    California oils and gases." Drilling and Production Practice. American
    Petroleum Institute, 1947.
    Uses refinment described here: https://petrowiki.org/Oil_bubblepoint_pressure
    :param density: density of oil at room conditions [kg/m^3]
    :param gas_oil_ratio: The volume ratio of gas to oil [l/l]
    :param gas_gravity: molar mass of gas relative to air molar mas.
    :param temperature: temperature of oil [°C]
    :return: bubble point pressure [Pa]
    """

    # Standing, M.B. (1947) uses:
    # * pressure in psi
    # * temperature in F
    # * gas oil ratio in cu ft per bbl
    # * bubble point in in absolute psi
    # * gravity_of_tank_oil in API
    #
    # The paper describes that bubble point is a function of
    # (gor / gr) ** 0.83
    #   * (10**(0.00091 * t)/10**(0.0125*gravity_of_tank_oil))
    #
    # and gives the following example:
    #   Bubble point pressure at 200° F
    #   of a liquid having gas-oil ratio
    #   350 CFB, a gas gravity 0.75, and
    #   a tank oil gravity of 30° API.
    #   The required pressure is found
    #   to be 1930 psia.
    #
    # We could scale to fit this example, however,
    # we use 1896 psia.
    #
    # For density in kg/m^3:
    # 10**(0.0125*gravity_of_tank_oil)=
    #        10**(1770.79/density - 1.64375)=
    #         e**(4072.69738323/density - 3.78487)
    # For temperature in celcius:
    #    10**(0.00091 * farenheit) =
    #    10**(0.001638*celcius + 0.02912) =
    #     e**(0.00377163*celsius + 0.0670513)
    #
    # Removing constants factors we get
    # (gor_ratio /
    #  gas_gravity)) ** 0.83 / ( np.exp(4072.69738323 / density -
    #  0.00377163438 * temperature_c)
    #
    # The equation occurs in this form as equation 21a in
    # Batzle, Michael, and Zhijing Wang. "Seismic properties of pore fluids."
    # Geophysics 57.11 (1992): 1396-1408.

    ratio = gas_oil_ratio / gas_gravity
    denominator = np.exp(4072 / density - 0.00377 * temperature)

    return 24469793.9134 * ratio ** 0.83 / denominator


def pressure_adjusted_dead_oil_density(pressure, reference_density):
    """
    Adjusts density of a dead oil (without dissolved gas) to a given pressure.

    Uses equation 18 from Batzle & Wang [1].

    :param reference_density: The density (g/cc) of the dead oil at 15.6 degrees celsius
        and atmospheric pressure.
    :param pressure: Pressure (MPa) to adjust to.
    :return: Density of oil at given pressure and 21 degrees celsius (~70 degrees farenheit).
    """
    return (
        reference_density
        + (0.00277 * pressure - 1.71 * 10 ** -7 * pressure ** 3)
        * (reference_density - 1.15) ** 2
        + 3.49 * 10 ** -4 * pressure
    )


def temperature_adjusted_dead_oil_density(temperature, density_at_21c):
    """
    Adjusts density of a dead oil (without dissolved gas) to a given temperature.

    Uses equation 19 from Batzle & Wang [1].

    :param density_at_21c: The density (g/cc) of the dead oil at 21 degrees celsius
    :param temperature: Temperature (celsius) of oil.
    :return: Density of oil at given temperature.
    """
    return density_at_21c / (0.972 + 3.81 * (10 ** -4) * (temperature + 17.78) ** 1.175)


def dead_oil_density(temperature, pressure, reference_density):
    """
    The density of oil without dissolved gas (dead).

    Uses equation 18 & 19 from Batzle & Wang [1].

    :param reference_density: Density of oil at 15.6 degrees celsius and atmospheric
        pressure. (g/cc)
    :param pressure: Pressure (MPa) of oil
    :param temperature: Temperature (celsius) of oil.
    :return: density of dead oil at given conditions.
    """
    density_p = pressure_adjusted_dead_oil_density(pressure, reference_density)
    return temperature_adjusted_dead_oil_density(temperature, density_p)


def dead_oil_primary_velocity(temperature, pressure, reference_density):
    """
    The primary wave velocity in oil without dissolved gas (dead).

    Uses equation 20a from Batzle & Wang [1].

    :param reference_density: Density of oil at 15.6 degrees celsius and atmospheric
        pressure (g/cc)
    :param pressure: Pressure (MPa) of oil
    :param temperature: Temperature (celsius) of oil.
    :return: primary velocity of dead oil in m/s.
    """
    return (
        2096 * np.sqrt(reference_density / (2.6 - reference_density))
        - 3.7 * temperature
        + 4.64 * pressure
        + 0.0115
        * (4.12 * np.sqrt(1.08 * reference_density ** -1 - 1) - 1)
        * temperature
        * pressure
    )


def dead_oil(temperature, pressure, reference_density):
    """
    :param reference_density: Density of the oil without dissolved gas
        at 15.6 degrees celsius and atmospheric pressure. kg/m3
    :param gas_oil_ratio: The volume ratio of gas to oil [l/l]
    :param gas_gravity: molar mass of gas relative to air molar mas.
    :param pressure: Pressure (Pa) of oil
    :param temperature: Temperature (celsius) of oil.
    :return: Material representing dead oil (without dissolved gas)
     as defined in Batzle & Wang [1].
    """
    return fluid(
        density=dead_oil_density(temperature, pressure * 1e-6, reference_density / 1000)
        * 1000,
        primary_velocity=dead_oil_primary_velocity(
            temperature, pressure * 1e-6, reference_density / 1000
        ),
    )


def live_oil_volume_factor(temperature, reference_density, gas_oil_ratio, gas_gravity):
    """
    Volume factor derived by Standing (1962), equation 23 in Batzle & Wang [1].
    :param reference_density: Density of the oil without dissolved gas
        at 15.6 degrees celsius and atmospheric pressure. (g/cc)
    :param gas_oil_ratio: The volume ratio of gas to oil [l/l]
    :param temperature: Temperature (celsius) of oil.
    :param gas_gravity: molar mass of gas relative to air molar mas.
    :return: A volume factor in calculating pseudo-density of live oil.
    """
    return (
        0.972
        + 0.00038
        * (
            2.4 * gas_oil_ratio * np.sqrt(gas_gravity / reference_density)
            + temperature
            + 17.8
        )
        ** 1.175
    )


def live_oil_pseudo_density(temperature, reference_density, gas_oil_ratio, gas_gravity):
    """
    Pseudo density used to substitute reference density in dead_oil_wave_velocity
    for live oils.

    Equation 22 in Batzle & Wang [1].

    :param reference_density: Density of the oil without dissolved gas
        at 15.6 degrees celsius and atmospheric pressure. (g/cc)
    :param pressure: Pressure (MPa) of oil
    :param gas_oil_ratio: The volume ratio of gas to oil [l/l]
    :param temperature: Temperature (celsius) of oil.
    :param gas_gravity: molar mass of gas relative to air molar mas.
    :return: Pseudo-density of live oil.
    """
    b0 = live_oil_volume_factor(
        temperature, reference_density, gas_oil_ratio, gas_gravity
    )
    return (reference_density / b0) / (1 + 0.001 * gas_oil_ratio)


def live_oil_density(
    temperature, pressure, reference_density, gas_oil_ratio, gas_gravity
):
    """
    Density of live oil at saturation.

    Equation 24 in Batzle & Wang [1].

    :param reference_density: Density of the oil without dissolved gas
        at 15.6 degrees celsius and atmospheric pressure. (g/cc)
    :param pressure: Pressure (MPa) of oil
    :param gas_oil_ratio: The volume ratio of gas to oil [l/l]
    :param temperature: Temperature (celsius) of oil.
    :param gas_gravity: molar mass of gas relative to air molar mas.
    :return: Density of live oil.
    """
    b0 = live_oil_volume_factor(
        temperature, reference_density, gas_oil_ratio, gas_gravity
    )
    return (reference_density + 0.0012 * gas_gravity * gas_oil_ratio) / b0


def live_oil_primary_velocity(
    temperature, pressure, reference_density, gas_oil_ratio, gas_gravity
):
    """
    Primary wave velocity of live oil at saturation.

    Substitute Equation 22 in Equation 20 of Batzle & Wang [1].

    :param reference_density: Density of the oil without dissolved gas
        at 15.6 degrees celsius and atmospheric pressure. (g/cc)
    :param pressure: Pressure (MPa) of oil
    :param gas_oil_ratio: The volume ratio of gas to oil [l/l]
    :param temperature: Temperature (celsius) of oil.
    :param gas_gravity: molar mass of gas relative to air molar mas.
    :return: Primary wave velocity of live oil.
    """
    rho_marked = live_oil_pseudo_density(
        temperature, reference_density, gas_oil_ratio, gas_gravity
    )
    return dead_oil_primary_velocity(temperature, pressure, rho_marked)


def live_oil(temperature, pressure, reference_density, gas_oil_ratio, gas_gravity):
    """
    :param reference_density: Density of the oil without dissolved gas
        at 15.6 degrees celsius and atmospheric pressure. (kg/m^3)
    :param gas_oil_ratio: The volume ratio of gas to oil [l/l]
    :param gas_gravity: molar mass of gas relative to air molar mas.
    :param pressure: Pressure (Pa) of oil
    :param temperature: Temperature (celsius) of oil.
    :return: Material representing live oil (with dissolved gas) as defined in Batzle & Wang [1].
    """
    if np.any(
        pressure
        < oil_bubble_point(reference_density, gas_oil_ratio, gas_gravity, temperature)
    ):
        warnings.warn(
            "Pressure is below bubble point of oil, estimated elastic properties can be inaccurate"
        )
    return fluid(
        density=live_oil_density(
            temperature,
            pressure * 1e-6,
            reference_density / 1000,
            gas_oil_ratio,
            gas_gravity,
        )
        * 1000,
        primary_velocity=live_oil_primary_velocity(
            temperature,
            pressure * 1e-6,
            reference_density / 1000,
            gas_oil_ratio,
            gas_gravity,
        ),
    )


def oil(temperature, pressure, reference_density, gas_oil_ratio, gas_gravity):
    """
    :param reference_density: Density of the oil without dissolved gas
        at 15.6 degrees celsius and atmospheric pressure. (kg/m^3)
    :param gas_oil_ratio: The volume ratio of gas to oil [l/l]
    :param gas_gravity: molar mass of gas relative to air molar mas.
    :param pressure: Pressure (Pa) of oil
    :param temperature: Temperature (celsius) of oil.
    :return: Material representing oil (with dissolved gas, or without if
        gas_oil_ratio=0.0) as defined in Batzle & Wang [1].
    """
    # Since live_oil with gas_oil_ratio=0.0 is not equal to dead oil
    # we use an apodization function to interpolate between the two

    def triangular_window(x, length=2):
        """
        A triangular window function around the origin, 1.0 at x=0.0, linear
        and 0.0 outside the window.
        :param length: total length of the window, ie., function is nonzero in
            [-length/2, length/2].
        :param x: which x to evaluate the window at
        :return: value of window function at x.
        """
        if x > length / 2:
            return 0.0
        if x < -length / 2:
            return 0.0
        return np.abs((np.abs(x) - length / 2) / (length / 2))

    loil = live_oil(
        temperature, pressure, reference_density, gas_oil_ratio, gas_gravity
    )
    doil = dead_oil(temperature, pressure, reference_density)
    window = triangular_window(gas_oil_ratio)
    return fluid(
        density=doil.density * window + (1 - window) * loil.density,
        primary_velocity=doil.primary_velocity * window
        + (1 - window) * loil.primary_velocity,
    )
