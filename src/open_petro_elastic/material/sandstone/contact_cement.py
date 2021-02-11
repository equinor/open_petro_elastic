"""
Dvorkin Nur (1996) contact cement model of cemented granular
material. Based on the assumption that porosity decreases from the
initial porosity of a sand pack due to the uniform
deposition of cement on the surface of the grains.

* Dvorkin, Jack, and Amos Nur. "Elasticity of high-porosity sandstones: Theory
for two North Sea data sets." Geophysics 61.5 (1996): 1363-1370.
"""

import numpy as np

from open_petro_elastic.material import Material


def shear_stiffness(
    cement_shear_modulus, grain_shear_modulus, grain_poisson_ratio, radius_ratio
):
    """
    Calculates the parameter St in (1) from dvorkin et. al (1991), also
    known as the tangential stiffness.

    See appendix A for description of calculation.

    :param grain_poisson_ratio: The poisson's ratio of the grains.
    :param cement_shear_modulus: The shear modulus of the cement.
    :param grain_shear_modulus: The shear modulus of the grain.
    :param radius_ratio: The ratio of the radius of the cement layer to the
        grain radius. Constant alpha in dvorkin et. al (1991) equations (1)-(3).
    """
    poiss = grain_poisson_ratio
    alpha = radius_ratio
    aat = cement_shear_modulus / (np.pi * grain_shear_modulus)
    at = (
        -1e-2
        * (2.26 * poiss ** 2 + 2.07 * poiss + 2.3)
        * aat ** (0.079 * poiss ** 2 + 0.1754 * poiss - 1.342)
    )
    bt = (0.0573 * poiss ** 2 + 0.0937 * poiss + 0.202) * aat ** (
        0.0274 * poiss ** 2 + 0.0529 * poiss - 0.8765
    )
    ct = (
        1e-4
        * (9.654 * poiss ** 2 + 4.945 * poiss + 3.1)
        * aat ** (0.01867 * poiss ** 2 + 0.4011 * poiss - 1.8186)
    )

    return at * alpha ** 2 + bt * alpha + ct


def normal_stiffness(
    grain_poisson_ratio,
    cement_poisson_ratio,
    cement_shear_modulus,
    grain_shear_modulus,
    radius_ratio,
):
    """
    Calculates the parameter Sn in (1) from dvorkin et. al (1991).
    See appendix A for description of calculation.
    :param grain_poisson_ratio: The poisson's ratio of the grains.
    :param cement_poisson_ratio: The poisson's ratio of the cement.
    :param cement_shear_modulus: The shear modulus of the cement.
    :param grain_shear_modulus: The shear modulus of the grain.
    :param radius_ratio: The ratio of the radius of the cement layer to the
        grain radius. Constant alpha in dvorkin et. al (1991) equations (1)-(3).
    """
    poiss = grain_poisson_ratio
    poiss_c = cement_poisson_ratio
    mu_cem = cement_shear_modulus
    mu_gr = grain_shear_modulus
    alpha = radius_ratio
    aan = (2 * mu_cem / (np.pi * mu_gr)) * (
        (1 - poiss) * (1 - poiss_c) / (1 - 2 * poiss_c)
    )
    cn = 0.00024649 * aan ** (-1.9864)
    bn = 0.20405 * aan ** (-0.89008)
    an = -0.024153 * aan ** (-1.3646)

    return an * alpha ** 2 + bn * alpha + cn


def cement_radius_ratio(cemented_sand_porosity, grain_porosity):
    """
    :return: the ratio of the radius of the cement layer to the grain
        radius, alpha in dvorkin et. al (1991) equation (2). Unlike the paper,
        the ratio is restricted to [0.0, 1.0] to avoid impossible values caused
        by numerical instability.
    :param cemented_sand_porosity: Porosity of cemented sand.
    :param grain_porosity: The porosity of the granular material.
    """
    if np.any(grain_porosity > 1) or np.any(grain_porosity < 0):
        raise ValueError(f"sand porosity must be between 0 and 1, was {grain_porosity}")
    if np.any(cemented_sand_porosity > 1) or np.any(cemented_sand_porosity < 0):
        raise ValueError(
            f"cemented sand porosity must be between 0 and 1, was {cemented_sand_porosity}"
        )
    if np.any(grain_porosity < cemented_sand_porosity):
        raise ValueError(
            f"sand porosity must be greater than cemented sand porosity: {grain_porosity} < {cemented_sand_porosity}"
        )
    calculated_value = (
        2 * (grain_porosity - cemented_sand_porosity) / (3 * (1 - grain_porosity))
    ) ** 0.5

    return np.clip(0.0, 1.0, calculated_value)


def contact_cement(
    cement,
    mineral,
    cemented_sand_porosity,
    sand_porosity=0.4,
    shear_reduction=1.0,
    coordination_number=9,
):
    """
    Cement is assumed to be evenly distributed on the
    grains (scheme 2 in Dvorkin and Nur's original paper)

    :param mineral: material representing the sand framework material.
    :param cement: material representing the pore-filling cement.
    :param cemented_sand_porosity: Porosity of cemented sand.
    :param sand_porosity: The porosity of the sand, defaults to 0.4.
        porosity of dry sand minus the porosity of cemented sand.
    :param coordination_number: Average number of contacts per grain,
        defaults to 9.
    :param shear_reduction: Shear reduction parameter used to account for
        tangential frictionaless grain contacts, defaults to no reduction, ie. 1.0.
    :return: A material representing the resulting cemented sand.
    """
    if np.any(sand_porosity > 1) or np.any(sand_porosity < 0):
        raise ValueError(f"sand porosity must be between 0 and 1, was {sand_porosity}")
    if np.any(cemented_sand_porosity > 1) or np.any(cemented_sand_porosity < 0):
        raise ValueError(
            f"cemented sand porosity must be between 0 and 1, was {cemented_sand_porosity}"
        )
    if np.any(sand_porosity < cemented_sand_porosity):
        raise ValueError(
            f"sand porosity must be greater than cemented sand porosity: {sand_porosity} < {cemented_sand_porosity}"
        )
    alpha = cement_radius_ratio(cemented_sand_porosity, sand_porosity)

    st = shear_stiffness(
        cement.shear_modulus, mineral.shear_modulus, mineral.poisson_ratio, alpha
    )
    sn = normal_stiffness(
        mineral.poisson_ratio,
        cement.poisson_ratio,
        cement.shear_modulus,
        mineral.shear_modulus,
        alpha,
    )

    # K_eff from equation (1)
    k_eff = (
        (1 / 6)
        * coordination_number
        * (1 - sand_porosity)
        * cement.primary_wave_modulus
        * sn
    )

    # G_eff from equations (1)
    mu_eff = (3 / 5) * k_eff + shear_reduction * (3 / 20) * coordination_number * (
        1 - sand_porosity
    ) * cement.shear_modulus * st

    return Material(
        bulk_modulus=k_eff,
        shear_modulus=mu_eff,
        density=(1 - cemented_sand_porosity)
        * (alpha * cement.density + (1 - alpha) * mineral.density),
    )
