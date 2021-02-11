import numpy as np

from ..hashin_shtrikman import hashin_shtrikman_walpole
from ..material import Material
from .constant_cement import constant_cement
from .contact_cement import contact_cement
from .friable_sand import friable_sand


def patchy_cement(
    sand,
    cement,
    porosity,
    contact_cement_porosity,
    upper_bound_porosity,
    pressure,
    lower_bound_pressure,
    critical_porosity=0.4,
    shear_reduction=1.0,
    coordination_number=9,
):
    """
    A hybrid sandstone model based on the patchy cement model of Avseth et. al (2016).

    :param sand: Material representing the sand of the sandstone.
    :param cement: Material representing the cement of the sandstone.
    :param porosity: Porosity of the uncemented sandstone.
    :param contact_cement_porosity: Porosity of the contact cement in the
        hybrid model.
    :param constant_cement_porosity: The porosity of constant cement
        in the hybrid model.
    :param pressure: Pressure used for friable sand.
    :param lower_bound_pressure: The pressure used for the lower bound
        friable sand.
    :param critical_porosity: The critical porosity of the sandstone.
    :param shear_reduction: Shear reduction factor.
    :param coordination number: The coordination number of the sand.

    Avseth, Per & Skjei, Norunn & Mavko, Gary. (2016). Rock-physics modeling of
    stress sensitivity and 4D time shifts in patchy cemented sandstones —
    Application to the Visund Field, North Sea. The Leading Edge. 35. 868-878.
    10.1190/tle35100868.1.
    """
    dense_packing = hashin_shtrikman_walpole(
        cement, sand, critical_porosity - upper_bound_porosity
    )

    friable = friable_sand(
        dense_packing,
        porosity,
        critical_porosity,
        pressure,
        coordination_number,
        shear_reduction,
    )
    lower_bound = friable_sand(
        dense_packing,
        porosity,
        critical_porosity,
        lower_bound_pressure,
        coordination_number,
        shear_reduction,
    )
    upper_bound = constant_cement(
        sand,
        cement,
        porosity,
        upper_bound_porosity,
        critical_porosity,
        coordination_number,
        shear_reduction,
    )
    contact = contact_cement(
        cement,
        sand,
        contact_cement_porosity,
        critical_porosity,
        shear_reduction,
        coordination_number,
    )
    constant = hashin_shtrikman_walpole(
        dense_packing, contact, 1 - porosity / contact_cement_porosity
    )

    kcc = np.asarray(constant.bulk_modulus)
    kf = np.asarray(friable.bulk_modulus)
    klow = np.asarray(lower_bound.bulk_modulus)
    kup = np.asarray(upper_bound.bulk_modulus)

    mucc = np.asarray(constant.shear_modulus)
    muf = np.asarray(friable.shear_modulus)
    mulow = np.asarray(lower_bound.shear_modulus)
    muup = np.asarray(upper_bound.shear_modulus)

    wk = np.ones(kcc.shape)
    wmu = np.ones(mucc.shape)

    idwk = kup == klow
    idwmu = muup == mulow

    wk[~idwk] = (kcc[~idwk] - klow[~idwk]) / (kup[~idwk] - klow[~idwk])
    wmu[~idwmu] = (mucc[~idwmu] - mulow[~idwmu]) / (muup[~idwmu] - mulow[~idwmu])

    wmu = np.clip(wmu, 0, 1)
    wk = np.clip(wk, 0, 1)

    cement_fraction = critical_porosity - contact_cement_porosity

    kps = kf + wk * (kup - kf)
    mups = muf + wmu * (muup - muf)
    rps = (
        1 - porosity - cement_fraction
    ) * sand.density + cement_fraction * cement.density

    return Material(bulk_modulus=kps, shear_modulus=mups, density=rps)
