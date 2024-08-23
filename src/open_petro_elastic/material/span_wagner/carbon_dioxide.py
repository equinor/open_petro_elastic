"""
Package for calculating CO2 properties based on Span & Wagner [2].
"""

import numpy as np
import pkg_resources
import scipy.optimize
from scipy.interpolate import RegularGridInterpolator

from open_petro_elastic import float_vectorize
from open_petro_elastic.material.fluid import fluid_material as fluid
from open_petro_elastic.material.span_wagner import equations
from open_petro_elastic.material.span_wagner.coefficients import a0, theta0
from open_petro_elastic.material.span_wagner.tables.lookup_table import (
    load_lookup_table_interpolator,
)

# Constants
CO2_GAS_CONSTANT = 0.1889241 * 1e3  # J / kg K
CO2_CRITICAL_TEMPERATURE = 304.1282  # K
CO2_CRITICAL_DENSITY = 467.6  # kg / m^3
CO2_CRITICAL_PRESSURE = 7.3773  # MPa
CO2_TRIPLE_TEMPERATURE = 216.592  # K
CO2_TRIPLE_PRESSURE = 0.51795  # MPa


def co2_helmholtz_energy(delta, tau, dd, dt):
    """
    Helmholtz energy as defined by equation 6.1 in Span & Wagner [2]

    :param delta: Reduced density, unit-less. That is, density / CO2_CRITICAL_DENSITY.
    :param tau: Inverse reduced temperature, unit-less. That is, CO2_CRITICAL_TEMPERATURE / (absolute) temperature
    :param dd: Degree of derivation wrt. delta. Integer between 0 and 2.
    :param dt: Degree of derivation wrt. tau. Integer between 0 and 2, as long as (dt + dd < 3)
    """
    return ideal_gas_helmholtz_energy(
        delta, tau, dd, dt
    ) + co2_residual_helmholtz_energy(delta, tau, dd, dt)


def ideal_gas_helmholtz_energy(delta, tau, dd, dt):
    """
    Helmholtz energy from ideal gas behavior as defined by equation 2.3 in Span & Wagner [2]. See function
    co2_helmholtz_energy for argument documentation.
    """
    # Adjust array shapes
    tau = np.asarray(tau)
    delta = np.asarray(delta)
    return_scalar = tau.ndim == 0
    tau2 = tau.reshape(-1, 1)  # Needed for array-based sums

    if dt == dd == 0:
        _sum = np.sum(a0[0, 3:] * np.log(1 - np.exp(-theta0 * tau2)), axis=-1)
        result = (
            np.log(delta) + a0[0, 0] + a0[0, 1] * tau + a0[0, 2] * np.log(tau) + _sum
        )
    elif dt == 1 and dd == 0:
        _sum = np.sum(
            a0[0, 3:] * theta0 * ((1 - np.exp(-theta0 * tau2)) ** -1 - 1), axis=-1
        )
        result = a0[0, 1] + a0[0, 2] / tau + _sum
    elif dt == 0 and dd == 1:
        return 1 / delta
    elif dt == 2 and dd == 0:
        _sum = np.sum(
            a0[0, 3:]
            * theta0**2
            * np.exp(-theta0 * tau2)
            * (1 - np.exp(-theta0 * tau2)) ** -2,
            axis=-1,
        )
        result = -a0[0, 2] / tau**2 - _sum
    elif dt == 0 and dd == 2:
        return -1 / delta**2
    elif dt == 1 and dd == 1:
        return 0
    else:
        raise NotImplementedError
    if return_scalar:
        return result[0]
    else:
        return result


def co2_residual_helmholtz_energy(delta, tau, dd, dt):
    """
    Residual part of Helmholtz energy as defined by the equation in Table 32 of Span & Wagner [2]. See
    co2_helmholtz_energy for argument documentation.
    """
    tau = np.asarray(tau)
    delta = np.asarray(delta)
    return_scalar = (tau.ndim == 0) & (delta.ndim == 0)
    tau = tau.reshape(-1, 1)
    delta = delta.reshape(-1, 1)
    # tau == 1.0 or delta == 1.0 leads to numerically invalid results. The values are nudged to avoid nan output.
    tau[tau == 1.0] -= 1e-15
    delta[delta == 1.0] -= 1e-15

    res = equations.residual_helmholtz_energy(delta, tau, dd, dt)

    if return_scalar:
        return res[0]
    else:
        return res


def carbon_dioxide_pressure(
    absolute_temperature, density, d_density=0, d_temperature=0, isentropic=False
):
    """
    CO2 pressure (MPa) as given by Table 3 of Span & Wagner [2]

    :param absolute_temperature: Temperature in K.
    :param density: CO2 density (kg / m^3).
    :param d_density: Degree of derivation wrt. density.
    :param d_temperature: Degree of derivation wrt. temperature.
    :param isentropic: Correction for isentropic conditions. Relevant primarily for isentropic bulk modulus
    """
    tau = CO2_CRITICAL_TEMPERATURE / absolute_temperature
    delta = density / CO2_CRITICAL_DENSITY
    if d_temperature != 0:
        raise NotImplementedError
    if d_density == 0:
        return (
            density
            * CO2_GAS_CONSTANT
            * absolute_temperature
            * (1 + delta * co2_residual_helmholtz_energy(delta, tau, 1, 0))
            / 1e6
        )
    elif d_density == 1:
        first = 2 * delta * co2_residual_helmholtz_energy(delta, tau, 1, 0)
        second = delta**2 * co2_residual_helmholtz_energy(delta, tau, 2, 0)
        if isentropic is False:
            third = 0
        else:
            # See Table 3 of Span & Wagner (speed of sound)
            nom = (
                1
                + delta * co2_residual_helmholtz_energy(delta, tau, 1, 0)
                - delta * tau * co2_residual_helmholtz_energy(delta, tau, 1, 1)
            ) ** 2
            den = tau**2 * (co2_helmholtz_energy(delta, tau, 0, 2))
            third = -nom / den
        return (
            absolute_temperature * CO2_GAS_CONSTANT * (1 + first + second + third) / 1e6
        )


def saturated_liquid_density(absolute_temperature):
    """
    Saturated liquid density as defined by equation 3.14 of Span & Wagner [2]

    :param absolute_temperature: Absolute temperature in K. Should satisfy:
        CO2_TRIPLE_TEMPERATURE < absolute_temperature < CO2_CRITICAL_TEMPERATURE
    """
    _a1 = 1.9245108
    _a2 = -0.62385555
    _a3 = -0.32731127
    _a4 = 0.39245142
    _t = 1 - absolute_temperature / CO2_CRITICAL_TEMPERATURE
    inner = _a1 * _t**0.34 + _a2 * _t**0.5 + _a3 * _t ** (10 / 6) + _a4 * _t ** (11 / 6)
    return CO2_CRITICAL_DENSITY * np.exp(inner)


def saturated_vapor_density(absolute_temperature):
    """
    Saturated vapor density as defined by equation 3.15 of Span & Wagner

    :param absolute_temperature: Absolute temperature in K. Should satisfy:
        CO2_TRIPLE_TEMPERATURE < absolute_temperature < CO2_CRITICAL_TEMPERATURE
    """
    # Assert temp < critical
    _a1 = -1.7074879
    _a2 = -0.82274670
    _a3 = -4.6008549
    _a4 = -10.111178
    _a5 = -29.742252
    _t = 1 - absolute_temperature / CO2_CRITICAL_TEMPERATURE
    inner = (
        _a1 * _t**0.34
        + _a2 * _t**0.5
        + _a3 * _t
        + _a4 * _t ** (7 / 3)
        + _a5 * _t ** (14 / 3)
    )
    return CO2_CRITICAL_DENSITY * np.exp(inner)


def sublimation_pressure(absolute_temperature):
    """
    Sublimation pressure as defined by equation 3.12 of Span & Wagner [2]

    :param absolute_temperature: Absolute temperature in K. Should satisfy absolute_temperature < CO2_TRIPLE_TEMPERATURE
    """
    _a1 = -14.740846
    _a2 = 2.4327015
    _a3 = -5.3061778
    _t = 1 - absolute_temperature / CO2_TRIPLE_TEMPERATURE
    inner = _a1 * _t + _a2 * _t**1.9 + _a3 * _t**2.9
    return CO2_TRIPLE_PRESSURE * np.exp(inner / (1 - _t))


def vapor_pressure(absolute_temperature):
    """
    Vapor pressure as defined by equation 3.13 of Span & Wagner [2]

    :param absolute_temperature: Absolute temperature in K. Should satisfy:
        CO2_TRIPLE_TEMPERATURE < absolute_temperature < CO2_CRITICAL_TEMPERATURE
    """
    _a1 = -7.0602087
    _a2 = 1.9391218
    _a3 = -1.6463597
    _a4 = -3.2995634
    _t = 1 - absolute_temperature / CO2_CRITICAL_TEMPERATURE
    inner = _a1 * _t**1.0 + _a2 * _t**1.5 + _a3 * _t**2.0 + _a4 * _t**4.0
    return CO2_CRITICAL_PRESSURE * np.exp(inner / (1 - _t))


def melting_pressure(absolute_temperature):
    """
    Melting pressure as defined by equation 3.10 of Span & Wagner [2]

    :param absolute_temperature: Absolute temperature in K. Should satisfy CO2_TRIPLE_TEMPERATURE < absolute_temperature
    """
    _a1 = 1955.5390
    _a2 = 2055.4593
    _t = absolute_temperature / CO2_TRIPLE_TEMPERATURE - 1
    return CO2_TRIPLE_PRESSURE * (1 + _a1 * _t + _a2 * _t**2)


def _determine_density_bounds(absolute_temperature, pressure, force_vapor):
    """
    Calculate the upper and lower bound on density
    """
    bounds = np.zeros((absolute_temperature.size, 2))
    bounds[:, 0] = 0.1
    bounds[:, 1] = 1500.0

    below_triple = absolute_temperature < CO2_TRIPLE_TEMPERATURE
    below_critical = ~below_triple & (absolute_temperature < CO2_CRITICAL_TEMPERATURE)

    bounds[below_triple, 1] = saturated_vapor_density(
        absolute_temperature[below_triple]
    )
    if force_vapor is True:
        bounds[below_critical, 1] = saturated_vapor_density(
            absolute_temperature[below_critical]
        )
    elif force_vapor is False:
        bounds[below_critical, 0] = saturated_liquid_density(
            absolute_temperature[below_critical]
        )
    else:  # force_vapor == 'auto'
        below_vapor_pressure = pressure < vapor_pressure(absolute_temperature)
        is_vapor = below_critical & below_vapor_pressure
        bounds[is_vapor, 1] = saturated_vapor_density(absolute_temperature[is_vapor])
        is_liquid = below_critical & ~below_vapor_pressure
        bounds[is_liquid, 0] = saturated_liquid_density(absolute_temperature[is_liquid])
    return bounds


def _find_initial_density_values(bounds, absolute_temperature, pressure):
    """
    Finds approximate density values for the provided temperature(s) and pressure(s). The result is only intended to be
    used by array_carbon_dioxide_density.
    """

    temps = np.geomspace(
        np.min(absolute_temperature) * 0.99, np.max(absolute_temperature) * 1.01, 41
    )
    press = np.geomspace(
        np.min(pressure) * 0.99, np.max(absolute_temperature) * 1.01, 41
    )
    tt, pp = np.meshgrid(temps, press, indexing="ij")
    densi = carbon_dioxide_density(
        tt.flatten(), pp.flatten(), force_vapor="auto", raise_error=False
    ).reshape(temps.size, press.size)
    rgi = RegularGridInterpolator(np.array((temps, press)), densi, method="linear")
    iv = rgi(np.array((absolute_temperature, pressure)).T)
    oob = (iv < bounds[:, 0]) | (iv > bounds[:, 1]) | np.isnan(iv)
    iv[oob] = np.mean(bounds[oob], axis=1)
    return iv


def array_carbon_dioxide_density(absolute_temperature, pressure, force_vapor):
    """
    Alternative implementation of a vectorized carbon dioxide density function. Implemented primarily for demonstration
    purposes. For large arrays, a look-up-table approach should be preferred.

    Utilizes scipy.optimize.newton, which is the only root-finding method of scipy that supports a vectorized functions.

    For argument documentation, see carbon_dioxide_density.
    """
    absolute_temperature = np.asarray(absolute_temperature)
    pressure = np.asarray(pressure)
    bounds = _determine_density_bounds(absolute_temperature, pressure, force_vapor)
    iv = _find_initial_density_values(bounds, absolute_temperature, pressure)
    opt = scipy.optimize.newton(
        lambda x: carbon_dioxide_pressure(absolute_temperature, x) - pressure,
        x0=iv,
        maxiter=10,
        full_output=True,
    )
    # scipy.optimize.newton may not always converge. We need to determine which of the elements of the solution are
    #  invalid. The opt.converged variable does not seem to suffice, so we perform separate checks. First, check that
    #  the solution is a valid root
    invalid = ~np.isclose(
        carbon_dioxide_pressure(absolute_temperature, opt.root),
        pressure,
        atol=1e-5,
        rtol=0.0,
    )

    # Next, check if the solution is anywhere out of bounds (since newton does not support brackets), and check for nan
    # values.
    invalid |= (
        (opt.root < bounds[:, 0]) | (opt.root > bounds[:, 1]) | np.isnan(opt.root)
    )

    # Finally, use the robust density method to determine the invalid results
    sol = opt.root
    sol[invalid] = carbon_dioxide_density(
        absolute_temperature[invalid], pressure[invalid], force_vapor=force_vapor
    )
    return sol


@float_vectorize
def _calculate_carbon_dioxide_density(
    absolute_temperature, pressure, force_vapor="auto", raise_error=True
):
    """
    Density of carbon dioxide. Found solving the Pressure equation of Table 3 in Span & Wagner [2] numerically for
    density. To ensure a single solution, the phase of the liquid must first be determined.

    :param absolute_temperature: Absolute temperature (K).
    :param pressure: Pressure (MPa).
    :param force_vapor: Boolean or 'auto'. If 'auto', the phase of the fluid is automatically determined. However, along
        the vaporization line (assuming T_triple < absolute_temperature < T_critical), the fluid is in two-phase
        equilibrium and the phase cannot be uniquely determined. If force_vapor is set to True, vapor phase is always
        selected, if False, liquid phase is selected. Outide the temperature bounds, this argument has no effect. This
        argument should only be used close to the vaporization boundary, otherwise the behavior might not be as
        expected.
    :param raise_error: Boolean. If True, raises an error if density cannot be determined. Otherwise, returns np.nan.

    :return: Density (kg / m^3)
    """
    bounds = _determine_density_bounds(
        np.array([absolute_temperature]), np.array([pressure]), force_vapor
    )[0, :]

    # Extend bounds slightly to ensure toms748 can converge
    bounds = [bounds[0] * 0.95, bounds[1] * 1.05]

    try:
        opt = scipy.optimize.root_scalar(
            lambda x: carbon_dioxide_pressure(absolute_temperature, x) - pressure,
            method="toms748",
            bracket=bounds,
            x0=np.sum(bounds) / 2,
        )
    except ValueError:
        if raise_error:
            raise
        else:
            return np.nan
    return opt.root


def carbon_dioxide_density(absolute_temperature, pressure, interpolate=False, **kwargs):
    """
    Density of carbon dioxide. Found either by direct calculation or interpolation. Any additional arguments are passed
    to _calculate_carbon_dioxide_density.

    :param absolute_temperature: Absolute temperature (K).
    :param pressure: Pressure (MPa).
    :param interpolate: Flag whether to interpolate data or not. If not, data is calculated directly. This is more
        accurate, but also more time-consuming. Data outside the bounds of the interpolator will be set to np.nan.

    :return: Density (kg / m^3)
    """
    if interpolate is False:
        return _calculate_carbon_dioxide_density(
            absolute_temperature, pressure, **kwargs
        )
    else:
        assert interpolate is True
        fp = pkg_resources.resource_filename(
            "open_petro_elastic.material.span_wagner.tables",
            "carbon_dioxide_density.npz",
        )
        interpolator = load_lookup_table_interpolator(fp)
        return interpolator(absolute_temperature, pressure)


def carbon_dioxide_bulk_modulus(absolute_temperature, density):
    """
    Isentropic bulk modulus, derived from the expression for speed of sound in Table 3 of Span & Wagner
    """
    d_pressure = carbon_dioxide_pressure(
        absolute_temperature, density, d_density=1, isentropic=True
    )
    return density * d_pressure * 1e6


def carbon_dioxide(absolute_temperature, pressure, density, **kwargs):
    """
    :param absolute_temperature: Temperature in Celsius. Valid range: 216 K - 1100 K
    :param pressure: Pressure in MPa. Valid range: 0 - 800
    :param density: Density in kg / m^3. May be provided instead of pressure. If both are provided, density has
        precedence
    """
    if density is None:
        density = carbon_dioxide_density(absolute_temperature, pressure, **kwargs)
    return fluid(
        density=density,
        bulk_modulus=carbon_dioxide_bulk_modulus(absolute_temperature, density),
    )
