"""
Models the elastic properties of a hydrostatic confining pressure applied to a
random indentical sphere packing, see Mavko et. al (2003) page 150

* "The Rock Physics Handbook: Tools for Seismic Analysis of Porous Media"
Gary Mavko, Tapan Mukerji, Jack Dvorkin
Cambridge University Press, 16. okt. 2003
"""
import numpy as np

from ..material import Material


def hertz_mindlin(
    mineral,
    porosity,
    pressure,
    shear_reduction=1.0,
    coordination_number=9,
):
    """
    Adjusts the moduli of a given granular material to a confining pressure
    according to Hertz-Mindlin model for pressure-induced moduli increase at
    critical porosity.

    :param mineral: The material composing the grains of a granular material.
    :param porisity: The porosity of the granular material.
    :param pressure: The confining pressure.
    :param shear_reduction: Shear reduction parameter used to account for
        tangential frictionaless grain contacts, defaults to no reduction, ie. 1.0.
    :param coordination_number: Average number of contacts per grain,
        defaults to 9.
    """
    v = mineral.poisson_ratio
    mu = mineral.shear_modulus
    n = coordination_number
    phi = porosity
    p = pressure

    a = ((3 * np.pi * (1 - v) * p) / (2 * n * (1 - phi) * mu)) ** (1 / 3)

    sn = (4 * a * mu) / (1 - v)
    st = (8 * a * mu) / (2 - v)

    k_eff = (n * (1 - phi) / (12 * np.pi)) * sn
    mu_eff = (n * (1 - phi) / (20 * np.pi)) * (sn + 1.5 * st * shear_reduction)

    return Material(
        bulk_modulus=k_eff,
        shear_modulus=mu_eff,
        density=mineral.density * (1 - porosity),
    )
