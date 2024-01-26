"""
Tests for open_petro_elastic.material.batzle_wang which referes to the
paper by Batzle & Wang (see help(open_petro_elastic.material.batzle_wang)).
Equations and figures in tests refer to that paper.

These figures are not completely accurate to computed values, possibly done by hand
and are sometimes from measured values.
"""

import pytest
from generators import positives, ratios
from hypothesis import assume, given
from predicates import assert_similar_material, material_between

import open_petro_elastic.material.batzle_wang as bw
from open_petro_elastic.material.batzle_wang.hydro_carbon_gas import (
    compressability_factor,
    compressability_rate_per_pseudoreduced_pressure,
)


@pytest.mark.parametrize(
    "gas_gravity, temperature, pressure, expected_density",
    [
        [0.6, 50, 10e6, 73],
        [0.6, 100, 10e6, 59],
        [0.6, 200, 10e6, 44],
        [0.6, 50, 25e6, 180],
        [0.6, 100, 25e6, 145],
        [0.6, 200, 25e6, 108],
        [0.6, 50, 50e6, 275],
        [0.6, 100, 50e6, 240],
        [0.6, 200, 50e6, 185],
        [1.2, 50, 10e6, 305],
        [1.2, 100, 10e6, 195],
        [1.2, 200, 10e6, 104],
        [1.2, 50, 25e6, 425],
        [1.2, 100, 25e6, 360],
        [1.2, 200, 25e6, 250],
        [1.2, 50, 50e6, 475],
        [1.2, 100, 50e6, 430],
        [1.2, 200, 50e6, 360],
    ],
)
def test_gas_density(gas_gravity, temperature, pressure, expected_density):
    """
    Tests that values calculated match that of figure 2 from Batzle & Wang.
    """
    assert bw.gas(temperature, pressure, gas_gravity).density == pytest.approx(
        expected_density, rel=0.01
    )


@pytest.mark.parametrize(
    "gas_gravity, temperature, pressure, expected_bulk_modulus",
    [
        [0.6, 50, 10e6, 175e5],
        [0.6, 100, 10e6, 182e5],
        [0.6, 200, 10e6, 19e6],
        [0.6, 50, 25e6, 57e6],
        [0.6, 100, 25e6, 525e5],
        [0.6, 200, 25e6, 515e5],
        [0.6, 50, 50e6, 163e6],
        [0.6, 100, 50e6, 129e6],
        [0.6, 200, 50e6, 104e6],
        [1.2, 50, 10e6, 303e5],
        [1.2, 100, 10e6, 18e6],
        [1.2, 200, 10e6, 175e5],
        [1.2, 50, 25e6, 207e6],
        [1.2, 100, 25e6, 125e6],
        [1.2, 200, 25e6, 62e6],
        [1.2, 50, 50e6, 576e6],
        [1.2, 100, 50e6, 340e6],
        [1.2, 200, 50e6, 182e6],
    ],
)
def test_gas_bulk_modulus(gas_gravity, temperature, pressure, expected_bulk_modulus):
    """
    Tests that values calculated match that of figure 3 from Batzle & Wang.
    """
    assert bw.gas(temperature, pressure, gas_gravity).bulk_modulus == pytest.approx(
        expected_bulk_modulus, rel=0.01
    )


@pytest.mark.parametrize(
    "reference_density, temperature, pressure, expected_density",
    [
        [780, 50, 0.1e6, 760],
        [780, 100, 0.1e6, 720],
        [780, 50, 25e6, 770],
        [780, 100, 25e6, 740],
        [780, 200, 25e6, 670],
        [780, 50, 50e6, 790],
        [780, 100, 50e6, 750],
        [780, 200, 50e6, 680],
        [880, 50, 0.1e6, 860],
        [880, 100, 0.1e6, 820],
        [880, 50, 25e6, 870],
        [880, 100, 25e6, 830],
        [880, 200, 25e6, 750],
        [880, 50, 50e6, 890],
        [880, 100, 50e6, 850],
        [880, 200, 50e6, 760],
        [1000, 50, 0.1e6, 970],
        [1000, 100, 0.1e6, 930],
        [1000, 50, 25e6, 980],
        [1000, 100, 25e6, 940],
        [1000, 200, 25e6, 850],
        [1000, 50, 50e6, 990],
        [1000, 100, 50e6, 950],
        [1000, 200, 50e6, 860],
    ],
)
def test_dead_oil_density(reference_density, temperature, pressure, expected_density):
    """
    Tests that values calculated match that of figure 5 from Batzle & Wang.
    """

    assert bw.dead_oil(
        temperature, pressure, reference_density
    ).density == pytest.approx(expected_density, rel=0.01)


@pytest.mark.parametrize(
    "reference_density, expected_velocity",
    [
        [1050, 1647.885],
        [950, 1513.2],
        [800, 1320.13],
        [700, 1195.04],
    ],
)
def test_dead_oil_velocity(reference_density, expected_velocity):
    """
    Tests that values calculated match that of figure 6 from Batzle & Wang.
    """

    atmospheric_pressure = 0.101325  # MPa
    room_temperature = 21.0
    assert bw.dead_oil(
        room_temperature, atmospheric_pressure, reference_density
    ).primary_velocity == pytest.approx(expected_velocity, rel=0.01)


@pytest.mark.parametrize(
    "salinity, temperature, expected_density",
    [
        [0.01e6, 0, 1007.47],
        [0.01e6, 10, 1007.07],
        [0.01e6, 20, 1005.34],
        [0.01e6, 25, 1004.09],
        [0.01e6, 30, 1002.61],
        [0.01e6, 40, 0999.08],
        [0.01e6, 50, 0994.82],
        [0.01e6, 60, 0990.00],
        [0.01e6, 80, 0978.50],
        [0.01e6, 100, 0965.10],
        [0.02e6, 0, 1015.09],
        [0.02e6, 10, 1014.42],
        [0.02e6, 20, 1012.46],
        [0.02e6, 25, 1011.12],
        [0.02e6, 30, 1009.57],
        [0.02e6, 40, 1005.93],
        [0.02e6, 50, 1001.61],
        [0.02e6, 60, 0996.70],
        [0.02e6, 80, 0985.20],
        [0.02e6, 100, 0971.90],
        [0.1e6, 0, 1076.77],
        [0.1e6, 10, 1074.19],
        [0.1e6, 20, 1070.68],
        [0.1e6, 25, 1068.79],
        [0.1e6, 30, 1066.76],
        [0.1e6, 40, 1062.38],
        [0.1e6, 50, 1057.53],
        [0.1e6, 60, 1052.30],
        [0.1e6, 80, 1040.50],
        [0.1e6, 100, 1027.60],
        [0.2e6, 0, 1156.63],
        [0.2e6, 10, 1152.54],
        [0.2e6, 20, 1147.79],
        [0.2e6, 25, 1145.33],
        [0.2e6, 30, 1142.85],
        [0.2e6, 40, 1137.74],
        [0.2e6, 50, 1132.38],
        [0.2e6, 60, 1126.80],
        [0.2e6, 80, 1114.60],
        [0.2e6, 100, 1101.70],
        [0.26e6, 0, 1207.09],
        [0.26e6, 10, 1202.54],
        [0.26e6, 20, 1197.17],
        [0.26e6, 25, 1194.43],
        [0.26e6, 30, 1191.70],
        [0.26e6, 40, 1186.14],
        [0.26e6, 50, 1180.45],
        [0.26e6, 60, 1174.70],
        [0.26e6, 80, 1162.60],
        [0.26e6, 100, 1149.20],
    ],
)
def test_brine_atmospheric_density_values(salinity, temperature, expected_density):
    """
    data points are from
    http://butane.chem.uiuc.edu/pshapley/genchem1/l21/1.html
    """
    atmospheric_pressure = 0.101325  # MPa
    assert bw.brine(
        temperature, atmospheric_pressure, salinity
    ).density == pytest.approx(expected_density, rel=0.006)


@pytest.mark.parametrize(
    "reference_density, temperature, pressure, expected_bulk_modulus",
    [
        [1000, 100, 0.1e6, 154e4],
        [1000, 100, 25e6, 186e4],
        [1000, 200, 25e6, 926e3],
        [1000, 100, 50e6, 221e4],
        [1000, 200, 50e6, 117e4],
        [880, 100, 0.1e6, 1045e3],
        [880, 100, 25e6, 1345e3],
        [880, 200, 25e6, 65e4],
        [880, 100, 50e6, 169e4],
        [880, 200, 50e6, 93e4],
        [780, 100, 0.1e6, 73e4],
        [780, 100, 25e6, 1e6],
        [780, 200, 25e6, 47e4],
        [780, 100, 50e6, 1325e3],
        [780, 200, 50e6, 75e4],
    ],
)
def test_dead_oil_bulk_modulus(
    reference_density, temperature, pressure, expected_bulk_modulus
):
    """
    Tests that values calculated match that of figure 7 from Batzle & Wang.
    """

    assert bw.dead_oil(
        temperature, pressure, reference_density
    ).bulk_modulus == pytest.approx(expected_bulk_modulus * 1000, rel=0.01)


@pytest.mark.parametrize(
    "temperature, pressure, expected_velocity",
    [
        [100, 50, 1647],
        [200, 50, 1480],
        [100, 100, 1734],
        [200, 100, 1620],
    ],
)
def test_water_velocity(temperature, pressure, expected_velocity):
    """
    Tests that values calculated match that of figure 12 from Batzle & Wang.

    This figure seem to be from data and not from computed values. Computed
    values seem to be imprecise for pressure above 100MPa.
    """

    assert bw.water(temperature, pressure).primary_velocity == pytest.approx(
        expected_velocity, rel=0.01
    )


@pytest.mark.parametrize(
    "salinity, temperature, pressure, expected_bulk_modulus",
    [
        [0.0, 50, 0.1e6, 2.35e9],
        [0.15e6, 50, 0.1e6, 3.14e9],
        [0.3e6, 50, 0.1e6, 4.13e9],
        [0.0, 100, 50e6, 2.66e9],
        [0.15e6, 100, 50e6, 3.39e9],
        [0.3e6, 100, 50e6, 4.26e9],
        [0.0, 300, 100e6, 1.54e9],
        [0.15e6, 300, 100e6, 2.16e9],
        [0.3e6, 300, 100e6, 3.12e9],
    ],
)
def test_brine_bulk_modulus(salinity, temperature, pressure, expected_bulk_modulus):
    """
    Tests that values calculated match that of figure 14 from Batzle & Wang.

    This figure seem to be from data and not from computed values. Computed
    values seem to be imprecise for pressure above 100MPa.
    """

    assert bw.brine(temperature, pressure, salinity).bulk_modulus == pytest.approx(
        expected_bulk_modulus, rel=0.01
    )


@pytest.mark.parametrize(
    "salinity, temperature, pressure, expected_density",
    [
        [0.0, 200, 9.81e6, 880],
        [0.0, 300, 49e6, 770],
        [0.02e6, 200, 49e6, 900],
        [0.02e6, 200, 98.1e6, 930],
        [0.15e6, 100, 9.81e6, 1070],
        [0.15e6, 100, 49e6, 1080],
        [0.15e6, 100, 98.1e6, 1090],
    ],
)
def test_brine_density(salinity, temperature, pressure, expected_density):
    """
    Tests that values calculated match that of figure 14 from Batzle & Wang.

    This figure seem to be from data and not from computed values. Computed
    values seem to be imprecise for pressure above 100MPa.
    """

    assert bw.brine(temperature, pressure, salinity).density == pytest.approx(
        expected_density, rel=0.01
    )


@pytest.mark.parametrize(
    "pressure, delta, epsilon",
    [(2, 0.01, 0.001), (2, 1, 0.01), (3, 1, 0.01), (4, 1, 0.01), (5, 2, 0.0001)],
)
def test_compressability_is_differentiated(pressure, delta, epsilon):
    """
    Tests that compressability_rate_per_pseudoreduced_pressure is
    compressability_factor differentiated. It is tested by

    f'(x) = lim_{h->0}{(f(x+h) - f(x)) / h}

    ie.

    for all epsilon, exists delta such that

    |((f(x+delta) - f(x)) / delta - f'(x)| < epsilon

    or (f(x+delta) - f(x)) / delta ~= f'(x)

    """

    def f(x):
        return compressability_factor(100, x, 0.6)

    def fd(x):
        """
        The differentiation is of pseudo reduced pressure so
        we use the chain rule here
        """
        return compressability_rate_per_pseudoreduced_pressure(100, x, 0.6) / (
            4.892 - 0.4048 * 0.6
        )

    def slope(x, delta):
        """
        Tangential slope of f at delta distance from x.
        """
        return (f(pressure + delta) - f(pressure)) / delta

    assert slope(pressure, delta) == pytest.approx(fd(pressure), abs=epsilon)


def below_bubble_point_gives_warning():
    with pytest.warns(UserWarning, match="Pressure is below bubble"):
        bw.oil(20.0, 1.0, 1.0, 200.0, 1.0)


pressures = positives(min_value=20e6, max_value=100e6)  # Pa
temperatures = positives(min_value=20, max_value=40)  # Celsius
reference_densities = positives(min_value=700, max_value=1000)  # kg/m3
gas_oil_ratios = positives(min_value=2.0, max_value=85.0)  # l/l
gas_gravities = positives(min_value=0.7, max_value=1.2)  # ratio


@given(pressures, temperatures, reference_densities, gas_gravities)
def test_oil_zero_gor_is_dead_oil(p, t, rd, g):
    assert_similar_material(bw.oil(t, p, rd, 0.0, g), bw.dead_oil(t, p, rd))


@given(pressures, temperatures, reference_densities, gas_oil_ratios, gas_gravities)
def test_oil_high_gor_is_live_oil(p, t, rd, gor, g):
    assert_similar_material(bw.oil(t, p, rd, gor, g), bw.live_oil(t, p, rd, gor, g))


@given(pressures, temperatures, reference_densities, ratios(), gas_gravities)
def test_oil_low_gor_is_between_live_and_dead(p, t, rd, gor, g):
    assert material_between(
        bw.live_oil(t, p, rd, gor, g), bw.dead_oil(t, p, rd), bw.oil(t, p, rd, gor, g)
    )


@given(pressures, temperatures, reference_densities, ratios(), ratios(), gas_gravities)
def test_oil_apodization_always_decreasing(p, t, rd, gor1, gor2, g):
    max_gor = max(gor1, gor2)
    min_gor = min(gor1, gor2)

    assume(max_gor - min_gor > 0.1)
    assume(bw.live_oil(t, p, rd, min_gor, g).density != bw.dead_oil(t, p, rd).density)

    min_diff = abs(
        bw.live_oil(t, p, rd, min_gor, g).density - bw.oil(t, p, rd, min_gor, g).density
    )
    max_diff = abs(
        bw.live_oil(t, p, rd, max_gor, g).density - bw.oil(t, p, rd, max_gor, g).density
    )

    assert max_diff < min_diff


@pytest.mark.parametrize(
    "gor, delta, epsilon",
    [
        (0.5, 1e-7, 10),
        (0.5, 1e-9, 0.1),
        (0.0, 1e-20, 1e-5),
        (0.9, 1e-15, 1e-3),
        (1.0, -1e-10, 0.01),
    ],
)
def test_oil_is_continous_on_unit_interval(gor, delta, epsilon):
    """
    Tests that oil is continous inside (0.0,1.0). It is

    tested by

    f(c) = lim_{x->c} f(x)

    ie.

    for all epsilon, exists delta such that

    |((f(c+delta) - f(c)| < epsilon

    or (f(c+delta) - f(c)) ~= f(c)

    """

    def f(x):
        return bw.oil(60, 10e6, 800, x, 0.6)

    assert_similar_material(f(gor), f(gor + delta), atol=epsilon, rtol=0.0)
