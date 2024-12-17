"""
Computation of bounds on elastic properties using the method
of Hashin-Shtrikman and their extensions.

see, for instance,

* Chinh, Pham Duc. "Bounds on the elastic moduli of statistically isotropic
  multicomponent materials and random cell polycrystals." International Journal
  of Solids and Structures 49.18 (2012): 2646-2659.
"""

from __future__ import annotations

import numpy as np

from open_petro_elastic.float_vectorize import float_vectorize

from .material import Material


def check_is_ratio(ratio: float, epsilon: float = 1e-5) -> None:
    """
    Check if the given ratio is less than 1.

    Parameters:
        ratio (float): The ratio to be checked.
        epsilon (float, optional): The tolerance value. Defaults to 1e-5.

    Raises:
        ValueError: If the ratio is greater than 1 + epsilon or less than 0 - epsilon.

    Returns:
        None
    """

    if np.any(1 + epsilon < ratio):
        raise ValueError(
            f"ratio in hashin-shtrikman bound should be less than 1: {ratio}"
        )
    if np.any(0 - epsilon > ratio):
        raise ValueError(
            f"ratio in hashin-shrtikman bound should be greater than 0: {ratio}"
        )


def hashin_shtrikman_bound(
    material1: Material, material2: Material, ratio: float
) -> Material:
    """
    Generate an instance of a Material object with Hashin-Shtrikman bounds for its bulk modulus, shear modulus, and density.

    Parameters:
        material1 (Material): The first constituent in the composite material.
        material2 (Material): The second constituent in the composite material.
        ratio (float): The volume fraction of material1 in the composite material.

    Returns:
        Material: The composite material with the calculated bulk modulus, shear modulus, and density.

    Raises:
        ValueError: If the ratio is not between 0 and 1.
        ValueError: If the material could not be created.

    References:
        - Hashin, Z., & Shtrikman, S. (1962).
        A variational approach to the theory of the elastic behavior of multiphase materials.
        Journal of the Mechanics and Physics of Solids, 10(4), 335-342.
    """
    check_is_ratio(ratio)
    k1 = material1.bulk_modulus
    k2 = material2.bulk_modulus
    mu1 = material1.shear_modulus
    mu2 = material2.shear_modulus
    f = ratio
    k = k1 + (1 - f) * (k2 - k1) / (1 + (k2 - k1) * f * (k1 + 4 / 3 * mu1) ** -1)
    mu = mu1 + (1 - f) * (mu2 - mu1) / (
        1 + 2 * (mu2 - mu1) * f * (k1 + 2 * mu1) / (5 * mu1 * (k1 + 4 / 3 * mu1))
    )

    return Material(
        bulk_modulus=k,
        shear_modulus=mu,
        density=ratio * material1.density + (1 - ratio) * material2.density,
    )


def hashin_shtrikman_average(
    material1: Material, material2: Material, ratio: float
) -> Material:
    """
    Generate an instance of a Material object with the Hashin-Shtrikman average properties.

    Parameters:
        material1 (Material): The first constituent in the composite material.
        material2 (Material): The second constituent in the composite material.
        ratio (float): The ratio of material1 to material2.

    Returns:
        Material: The material with the Hashin-Shtrikman average properties.

    Raises:
        ValueError: If the ratio is not between 0 and 1,
        ValueError: If there are errors in creating the bound materials.
        ValueError: If the average material instance cannot be created.

    Examples:
        >>> material1 = Material(bulk_modulus=100, shear_modulus=50, density=2.5)
        >>> material2 = Material(bulk_modulus=200, shear_modulus=110, density=3.0)
        >>> ratio = 0.6
        >>> result = hashin_shtrikman_average(material1, material2, ratio)
        >>> result.bulk_modulus
        130.79283887468029
        >>> result.shear_modulus
        68.78970958579256
        >>> result.density
        2.7
    """
    check_is_ratio(ratio)

    bound1 = hashin_shtrikman_bound(material1, material2, ratio)
    bound2 = hashin_shtrikman_bound(material2, material1, 1 - ratio)

    return Material(
        bulk_modulus=(bound1.bulk_modulus + bound2.bulk_modulus) / 2,
        shear_modulus=(bound1.shear_modulus + bound2.shear_modulus) / 2,
        density=ratio * material1.density + (1 - ratio) * material2.density,
    )


def hashin_shtrikman_walpole(
    material1: Material, material2: Material, ratio: float, bound: str = "lower"
) -> Material:
    """
    Generate an instance of a material at the lower(upper) bound corresponding to the
    Walpole refinement of the Hashin-Shtrikman bounds.

    Parameters:
        material1 (Material): The first constituent in the composite material.
        material2 (Material): The second constituent in the composite material.
        ratio (float): The ratio of material1 in the composite material.
        bound (str, optional): Either "lower" or "upper" bound. Defaults to "lower".

    Returns:
        Material: A material with lower (or upper) Hashin-Shtrikman-Walpole bound.

    Raises:
        ValueError: If the ratio is not between 0 and 1.
        ValueError: If the material could not be created.

    References:
        - Chinh, Pham Duc.
        "Bounds on the elastic moduli of statistically isotropic multicomponent materials and random cell polycrystals."
        International Journal of Solids and Structures 49.18 (2012): 2646-2659.
    """
    check_is_ratio(ratio)
    if bound not in ["lower", "upper"]:
        raise ValueError(f"bound should be either 'lower' or 'upper', got {bound}")

    @float_vectorize
    def calc_density(r1, r2, f1, f2):
        if r1 == r2:
            # Densities are equal, combined material has same density as either
            return r2
        return f1 * r1 + f2 * r2

    @float_vectorize
    def calc_bulk_modulus(k1, k2, mum, f1, f2):
        if k1 == k2:
            # Bulk moduli are equal, combined material has same modulus as either
            return k1

        return k1 + f2 / ((k2 - k1) ** -1 + f1 * (k1 + 4 / 3 * mum) ** -1)

    @float_vectorize
    def calc_shear_modulus(mu1, mu2, km, mum, f1, f2):
        if mu1 == mu2:
            # Shear moduli are equal, combined material has same modulus as either
            return mu1
        return mu1 + f2 / (
            (mu2 - mu1) ** -1
            + f1 * (mu1 + mum / 6 * (9 * km + 8 * mum) / (km + 2 * mum)) ** -1
        )

    ratio2 = 1 - ratio
    if bound == "lower":
        km = np.minimum(material1.bulk_modulus, material2.bulk_modulus)
        mum = np.minimum(material1.shear_modulus, material2.shear_modulus)
    else:
        km = np.maximum(material1.bulk_modulus, material2.bulk_modulus)
        mum = np.maximum(material1.shear_modulus, material2.shear_modulus)

    return Material(
        density=calc_density(material1.density, material2.density, ratio, ratio2),
        bulk_modulus=calc_bulk_modulus(
            material1.bulk_modulus, material2.bulk_modulus, mum, ratio, ratio2
        ),
        shear_modulus=calc_shear_modulus(
            material1.shear_modulus, material2.shear_modulus, km, mum, ratio, ratio2
        ),
    )
