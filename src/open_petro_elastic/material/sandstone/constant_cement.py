"""
The constant cement model (Aveset et. al (2000)) which assumes that sands of
varying porosity all have the same amount of contact cement.

* Avseth, Per, et al. "Rock physics diagnostic of North Sea sands: Link between
microstructure and seismic properties." Geophysical Research Letters 27.17
(2000): 2761-2764.
"""

from ..hashin_shtrikman import hashin_shtrikman_walpole
from .contact_cement import contact_cement


def constant_cement(
    sand,
    cement,
    porosity,
    contact_cement_porosity,
    sand_porosity=0.4,
    coordination_number=9,
    shear_reduction=1.0,
):
    """
    :param sand: material representing the sand framework material.
    :param cement: material representing a pore-filling cement.
    :param porosity: The porosity of constant cement.
    :param contact_cement_porosity: Porosity of the contact cement.
    :param sand_porosity: The porosity of the sand, defaults to 0.4.
        porosity of dry sand minus the porosity of cemented sand.
    :param coordination_number: Average number of contacts per grain,
        defaults to 9.
    :param shear_reduction: Shear reduction parameter used to account for
        tangential frictionaless grain contacts, defaults to no reduction, ie. 1.0.
    """

    # At the zero-porosity point, all original porosity is filled with grains.
    # The cement fraction surrounds the original grains.
    dense_packing = hashin_shtrikman_walpole(
        cement, sand, sand_porosity - contact_cement_porosity
    )

    # Dry rock properties of high-porosity end member calculated with
    # Dvorkin-Nur equation.
    cc = contact_cement(
        cement,
        sand,
        contact_cement_porosity,
        sand_porosity,
        shear_reduction,
        coordination_number,
    )

    # Fraction of zero-porosity end member
    fraction = 1 - porosity / contact_cement_porosity
    return hashin_shtrikman_walpole(dense_packing, cc, fraction)
