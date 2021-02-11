from ..hashin_shtrikman import hashin_shtrikman_walpole
from .hertz_mindlin import hertz_mindlin


def friable_sand(
    mineral,
    porosity,
    critical_porosity,
    effective_pressure,
    coordination_number=9,
    shear_reduction=1.0,
):
    """
    Models sandstone as a hashin-shtrikman-walpole lower bound
    of hertz_mindlin granular material at critical porosity and
    solid mineral.

    Described in Dvorkin, Jack, and Amos Nur. "Elasticity of high-porosity
    sandstones: Theory for two North Sea data sets." Geophysics 61.5 (1996): 1363-1370.

    :param mineral: The (composite) mineral the sandstone is made up of.
    :param porosity: The porosity of the sandstone.
    :param critical_porosity: The critical_porosity of the mineral.
    :param effective_pressure: The pressure used for hertz-mindlin grain.
    :param shear_reduction: Shear reduction parameter used to account for
        tangential frictionaless grain contacts, defaults to no reduction, ie. 1.0.
    :param coordination_number: Average number of contacts per grain,
        defaults to 9.
    """

    # Dry rock properties of high-porosity end member calculated with
    # Hertz-Mindlin equation
    critical_mineral = hertz_mindlin(
        mineral,
        critical_porosity,
        effective_pressure,
        shear_reduction,
        coordination_number,
    )

    # Fraction of solid
    fraction = 1 - (porosity / critical_porosity)

    # Hashin-Shtrikman lower bound describes the dry rock property mixing from
    # mineral properties to high-end porosity.
    return hashin_shtrikman_walpole(mineral, critical_mineral, fraction)
