"""
Recreates the figures from Batzle & Wang (see help(open_petro_elastic.material.batzle_wang)).
Equations and figures in tests refer to that paper.

These figures are not completely accurate to computed values, possibly done by hand
and are sometimes from measured values.
"""

import matplotlib.pyplot as plt
import numpy as np

import open_petro_elastic.material.batzle_wang as bw


def write_figure2():
    """
    Creates figure 2 from Batzle & Wang using the gas_density function, and
    writes it to figure_2.png.
    """
    fig, ax = plt.subplots()
    ax.set_ylim(0.0, 0.6)

    ax.set_xlabel("TEMPERATURE (°C)")
    ax.set_ylabel("GAS DENSITY (g/cm³)")
    xs = np.linspace(0, 350, 1000)

    def plot_gas_density(gas_gravity, pressure, *args, **kwargs):
        return ax.plot(
            xs,
            [bw.gas(x, pressure * 1e6, gas_gravity).density / 1e3 for x in xs],
            *args,
            **kwargs
        )

    plot_gas_density(0.6, 10, "k--", label="G = 0.6")
    plot_gas_density(1.2, 10, "k-", label="G = 1.2")
    plot_gas_density(1.2, 25, "k-")
    plot_gas_density(1.2, 50, "k-")
    plot_gas_density(0.6, 25, "k--")
    plot_gas_density(0.6, 50, "k--")

    ax.text(15.0, 0.09, "10", rotation=-20)
    ax.text(15.0, 0.21, "25", rotation=-32)
    ax.text(15.0, 0.265, "50", rotation=-21)

    ax.text(105.0, 0.18, "10", rotation=-31)
    ax.text(100.0, 0.35, "25", rotation=-32)
    ax.text(100.0, 0.415, "50MPa", rotation=-21)

    ax.legend()

    fig.savefig("figure_2.png")


def write_figure3():
    """
    Creates figure 3 from Batzle & Wang using the gas_bulk_modulus function
    and writes it to figure_3.png.
    """
    fig, ax = plt.subplots()
    ax.set_ylim(0.0, 650.0)
    ax.set_xlabel("TEMPERATURE (°C)")
    ax.set_ylabel("GAS BULK MODULUS (MPa)")
    xs = np.linspace(0.0, 350, 1000)
    xs_short = np.linspace(40, 350, 1000)

    def plot_gas_bulk_modulus(gas_gravity, temperatures, pressure, *args, **kwargs):
        return ax.plot(
            temperatures,
            [
                bw.gas(x, pressure * 1e6, gas_gravity).bulk_modulus / 1e6
                for x in temperatures
            ],
            *args,
            **kwargs
        )

    plot_gas_bulk_modulus(0.6, xs, 10, "k--", label="G = 0.6")
    plot_gas_bulk_modulus(1.2, xs, 10, "k-", label="G = 1.2")
    plot_gas_bulk_modulus(1.2, xs, 25, "k-")
    plot_gas_bulk_modulus(1.2, xs_short, 50, "k-")
    plot_gas_bulk_modulus(0.6, xs, 25, "k--")
    plot_gas_bulk_modulus(0.6, xs, 50, "k--")

    ax.text(100.0, 25, "10", rotation=0)
    ax.text(70.0, 60, "25", rotation=-2)
    ax.text(45.0, 120, "50", rotation=-18)

    ax.text(8.0, 22, "10", rotation=-31)
    ax.text(22.0, 300, "25", rotation=-70)
    ax.text(50.0, 500, "50MPa", rotation=-70)

    ax.legend()

    fig.savefig("figure_3.png")


def write_figure5():
    """
    Creates figure 5 from Batzle & Wang using the dead_oil_density function
    and writes it to figure_5.png
    """
    fig, ax = plt.subplots()
    ax.set_ylim(0.55, 1.05)
    ax.set_xlabel("TEMPERATURE (°C)")
    ax.set_ylabel("OIL DENSITY (g/cm³)")
    xs = np.linspace(0.0, 350, 1000)
    xs_short = np.linspace(0.0, 100, 1000)

    def plot_oil_density(reference_density, temperatures, pressure, *args, **kwargs):
        return ax.plot(
            temperatures,
            [
                bw.dead_oil(x, pressure * 1e6, reference_density * 1000).density / 1000
                for x in temperatures
            ],
            *args,
            **kwargs
        )

    plot_oil_density(1.00, xs_short, 0.1, "k-.", label="10 deg. API")
    plot_oil_density(0.88, xs_short, 0.1, "k-", label="30 deg. API")
    plot_oil_density(0.78, xs_short, 0.1, "k--", label="50 deg. API")

    plot_oil_density(1.00, xs, 25, "k-.")
    plot_oil_density(0.88, xs, 25, "k-")
    plot_oil_density(0.78, xs, 25, "k--")

    plot_oil_density(1.00, xs, 50, "k-.")
    plot_oil_density(0.88, xs, 50, "k-")
    plot_oil_density(0.78, xs, 50, "k--")

    ax.text(38.0, 0.72, "0.1 MPa", rotation=-18)
    ax.text(150.0, 0.66, "25 MPa", rotation=-18)
    ax.text(220.0, 0.66, "50 MPa", rotation=-18)

    ax.legend()

    fig.savefig("figure_5.png")


def write_figure6():
    """
    Creates figure 6 from Batzle & Wang using the dead_oil_wave_velocity function
    and writes it to figure_6.png.
    """
    fig, ax = plt.subplots()
    ax.set_ylim(1.1, 1.8)
    ax.set_xlim(1.1, 0.69)
    ax.set_xlabel("ρ₀")
    ax.set_ylabel("VELOCITY (km/sec)")
    xs = np.linspace(1.1, 0.7, 1000)
    ax.plot(
        xs,
        [bw.dead_oil(15.1, 0.1e6, x * 1e3).primary_velocity / 1000 for x in xs],
        "k-",
    )

    fig.savefig("figure_6.png")


def write_figure7():
    """
    Creates figure 7 from Batzle & Wang using the dead_oil density and velocity
    functions and writes it to figure_7.png.
    """
    fig, ax = plt.subplots()
    ax.set_xlabel("TEMPERATURE (°C)")
    ax.set_ylabel("OIL BULK MODULUS (MPa)")
    ax.set_ylim(0.0, 3400)
    xs = np.linspace(0.0, 350, 1000)
    xs_short = np.linspace(0.0, 100, 1000)

    def plot_oil_bulk_modulus(
        reference_density, temperatures, pressure, *args, **kwargs
    ):
        return ax.plot(
            temperatures,
            [
                bw.dead_oil(x, pressure * 1e6, reference_density * 1000).bulk_modulus
                / 1e6
                for x in temperatures
            ],
            *args,
            **kwargs
        )

    plot_oil_bulk_modulus(
        1.00,
        xs_short,
        0.1,
        "k-.",
        label="10 deg. API",
    )
    plot_oil_bulk_modulus(
        0.88,
        xs_short,
        0.1,
        "k-",
        label="30 deg. API",
    )
    plot_oil_bulk_modulus(
        0.78,
        xs_short,
        0.1,
        "k--",
        label="50 deg. API",
    )
    plot_oil_bulk_modulus(0.78, xs, 25, "k--")
    plot_oil_bulk_modulus(0.78, xs, 50, "k--")
    plot_oil_bulk_modulus(0.88, xs, 25, "k-")
    plot_oil_bulk_modulus(0.88, xs, 50, "k-")
    plot_oil_bulk_modulus(1.0, xs, 25, "k-.")
    plot_oil_bulk_modulus(1.0, xs, 50, "k-.")

    ax.legend()

    fig.savefig("figure_7.png")


def write_figure12():
    """
    Writes figure 12 from Batzle & Wang using the brine density and velocity functions
    and writes it to figure_12.png
    """
    fig, ax = plt.subplots()
    ax.set_ylim(0.5, 2.0)
    ax.set_xlim(0, 350)
    ax.set_xlabel("TEMPERATURE (°C)")
    ax.set_ylabel("VELOCITY (km/s)")
    xs = np.linspace(0.0, 350, 1000)

    ax.plot(
        xs,
        [bw.water(x, 1.0).primary_velocity / 1000 for x in xs],
        "k-",
        xs,
        [bw.water(x, 50).primary_velocity / 1000 for x in xs],
        "k-",
        xs,
        [bw.water(x, 100).primary_velocity / 1000 for x in xs],
        "k-",
    )

    fig.savefig("figure_12.png")


def write_figure8():
    """
    Creates figure 8 from Batzle & Wang using the dead_oil density and velocity
    functions and writes it to figure_8.png.
    """
    fig, ax = plt.subplots()
    ax.set_xlabel("PRESSURE (°C)")
    ax.set_ylabel("VELOCITY (km/s)")
    ax.set_ylim(0.9, 1.7)
    ax.set_xlim(0.0, 40.0)
    xs = np.linspace(1.0, 40, 1000)
    xs_short = np.linspace(5.0, 40, 1000)

    def plot_dead_oil_primary_velocity(temperature, pressures, *args, **kwargs):
        return ax.plot(
            pressures,
            [
                bw.dead_oil(temperature, x * 1e6, 916).primary_velocity / 1e3
                for x in pressures
            ],
            *args,
            **kwargs
        )

    def plot_oil_primary_velocity(temperature, pressures, gor, *args, **kwargs):
        return ax.plot(
            pressures,
            [
                bw.oil(temperature, x * 1e6, 916, gor, 0.8).primary_velocity / 1e3
                for x in pressures
            ],
            *args,
            **kwargs
        )

    plot_dead_oil_primary_velocity(
        22.8,
        xs,
        "k--",
    )
    plot_dead_oil_primary_velocity(
        72.0,
        xs,
        "k--",
    )

    plot_oil_primary_velocity(
        22.8,
        xs_short,
        85.0,
        "k--",
    )
    plot_oil_primary_velocity(
        72.0,
        xs_short,
        85.0,
        "k--",
    )

    ax.legend()

    fig.savefig("figure_8.png")


def write_figure13():
    """
    Writes figure 13 from Batzle & Wang using the brine density and velocity functions
    and writes it to figure_13.png
    """
    fig, ax = plt.subplots()
    ax.set_ylim(0.7, 1.3)
    ax.set_xlim(20, 351)
    ax.set_xlabel("TEMPERATURE (°C)")
    ax.set_ylabel("DENSITY (g/cc)")
    xs = np.linspace(0.0, 350, 1000)

    def plot_brine_density(salinity, pressure, *args, **kwargs):
        return ax.plot(
            xs,
            [bw.brine(x, pressure * 1e6, salinity).density / 1000 for x in xs],
            *args,
            **kwargs
        )

    plot_brine_density(0.0, 9.81, "k--")
    plot_brine_density(0.0, 49, "k--")
    plot_brine_density(0.0, 98.1, "k--")

    for salinity in [0.02e6, 0.15e6, 0.24e6]:
        plot_brine_density(salinity, 9.81, "k-")
        plot_brine_density(salinity, 49, "k-")
        plot_brine_density(salinity, 98.1, "k-")

    fig.savefig("figure_13.png")


def write_figure14():
    """
    Writes figure 14 from Batzle & Wang using the brine density and velocity functions
    and writes it to figure_14.png
    """
    fig, ax = plt.subplots()
    ax.set_xlim(20, 350)
    ax.set_xlabel("TEMPERATURE (°C)")
    ax.set_ylabel("BULK MODULUS (GPa)")
    xs = np.linspace(0.0, 350, 1000)
    xs_short = np.linspace(0.0, 100, 1000)

    def plot_brine_bulk_modulus(salinity, temperatures, pressure, *args, **kwargs):
        return ax.plot(
            temperatures,
            [
                bw.brine(x, pressure * 1e6, salinity * 1e6).bulk_modulus / 1e9
                for x in temperatures
            ],
            *args,
            **kwargs
        )

    plot_brine_bulk_modulus(
        0.0,
        xs_short,
        0.1,
        "k--",
        label="PPM = 0",
    )
    plot_brine_bulk_modulus(
        0.15,
        xs_short,
        0.1,
        "k-",
        label="PPM = 150000",
    )
    plot_brine_bulk_modulus(
        0.3,
        xs_short,
        0.1,
        "k-.",
        label="PPM = 300000",
    )

    plot_brine_bulk_modulus(0.0, xs, 50, "k--")
    plot_brine_bulk_modulus(0.15, xs, 50, "k-")
    plot_brine_bulk_modulus(0.3, xs, 50, "k-.")
    plot_brine_bulk_modulus(0.0, xs, 100, "k--")
    plot_brine_bulk_modulus(0.15, xs, 100, "k-")
    plot_brine_bulk_modulus(0.3, xs, 100, "k-.")

    ax.legend()

    fig.savefig("figure_14.png")


def main():
    write_figure2()
    write_figure3()
    write_figure5()
    write_figure6()
    write_figure7()
    write_figure8()
    write_figure12()
    write_figure13()
    write_figure14()


if __name__ == "__main__":
    main()
