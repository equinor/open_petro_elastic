"""
Computation of bounds on elastic properties using the method
of Hashin-Shtrikman and their extensions.

see, for instance,
* Chinh, Pham Duc. "Bounds on the elastic moduli of statistically isotropic
    multicomponent materials and random cell polycrystals." International Journal
    of Solids and Structures 49.18 (2012): 2646-2659.
"""
import numpy as np

from open_petro_elastic.float_vectorize import float_vectorize

from .material import Material


def check_is_ratio(ratio, epsilon=1e-5):
    """
    Checks that -epsilon <= ratio <= 1.0 + epsilon, throws ValueError if not.
    """
    if np.any(1 + epsilon < ratio):
        raise ValueError(
            f"ratio in hashin-shtrikman bound should be less than 1: {ratio}"
        )
    if np.any(0 - epsilon > ratio):
        raise ValueError(
            f"ratio in hashin-shrtikman bound should be greater than 0: {ratio}"
        )


def hashin_shtrikman_bound(material1, material2, ratio):
    """
    Gives a bound on moduli for composite materials.

    see, for instance, https://www.subsurfwiki.org/wiki/Hashin-Shtrikman_bounds
    :param material1: The first constiuent in the composite material.
    :param material2: The second constiuent in the composite material.
    :param ratio: The ratio of material1 in the composite material.
    :returns: The hashin-shrikman bound on the composite material. Its a
        upper bound if material1.bulk_modulus > material2.bulk_modulus and
        lower bound otherwise.
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


def hashin_shtrikman_average(material1, material2, ratio):
    """
    Average of Hashin-Shtrikman upper and lower bound.

    :param material1: The first constiuent in the composite material.
    :param material2: The second constiuent in the composite material.
    :param ratio: The ratio of material1 in the composite material.
    :returns: The average of upper and lower hashin_shtrikman bounds.
    """
    check_is_ratio(ratio)
    bound1 = hashin_shtrikman_bound(material1, material2, ratio)
    bound2 = hashin_shtrikman_bound(material2, material1, 1 - ratio)

    return Material(
        bulk_modulus=(bound1.bulk_modulus + bound2.bulk_modulus) / 2,
        shear_modulus=(bound1.shear_modulus + bound2.shear_modulus) / 2,
        density=ratio * material1.density + (1 - ratio) * material2.density,
    )


def hashin_shtrikman_walpole(material1, material2, ratio, bound="lower"):
    """
    Refinement of hashin_shtrikman bounds
    Chinh, Pham Duc. "Bounds on the elastic moduli of statistically isotropic
        multicomponent materials and random cell polycrystals." International
        Journal of Solids and Structures 49.18 (2012): 2646-2659.
    :param material1: The first constiuent in the composite material.
    :param material2: The second constiuent in the composite material.
    :param ratio: The ratio of material1 in the composite material.
    :param bound: Either "lower" or "upper", whether to return the
        lower or upper bound.
    :return: Lower (or upper) hashin-shtrikman-walpole bound for mixing
        of multicomponent materials.
    """
    check_is_ratio(ratio)

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
        mum = np.minimum(material2.shear_modulus, material2.shear_modulus)
    else:
        km = np.maximum(material1.bulk_modulus, material2.bulk_modulus)
        mum = np.maximum(material2.shear_modulus, material2.shear_modulus)

    return Material(
        density=calc_density(material1.density, material2.density, ratio, ratio2),
        bulk_modulus=calc_bulk_modulus(
            material1.bulk_modulus, material2.bulk_modulus, mum, ratio, ratio2
        ),
        shear_modulus=calc_shear_modulus(
            material1.shear_modulus, material2.shear_modulus, km, mum, ratio, ratio2
        ),
    )
