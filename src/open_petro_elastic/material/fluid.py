"""
Fluid is short hand for a material with shear modulus equal to zero.
"""

import numpy as np

from open_petro_elastic.float_vectorize import float_vectorize

from .material import Material


def fluid_material(*args, **kwargs):
    """
    fluid is a material with shear modulus=0
    """
    return Material(shear_modulus=0, *args, **kwargs)  # noqa: B026


def sum_saturations_to_one(saturations):
    total_saturation = sum(saturations)
    adjusted_saturation = np.maximum(total_saturation, 1e-15)
    return [
        np.where(total_saturation > 0.0, s, 1e-15) / adjusted_saturation
        for s in saturations
    ]


def mix_densities(fluids, saturations):
    return sum(
        fluid.density * saturation for (fluid, saturation) in zip(fluids, saturations)
    )


def mix_bulk_moduli(fluids, saturations):
    return 1.0 / sum(s / f.bulk_modulus for s, f in zip(saturations, fluids))


def wood_fluid_mixing(fluids, saturations):
    """
    woods function for fluid mixing, also known as Reuss average.
    :param fluids: A list of liquids to mix with gas.
    :param saturations: List of saturation for each liquid (sums to 1.0).
    :return: Material representing the mixed fluids.

    * A. B. Wood, A Textbook of Sound, 1st ed. MacMillan, New York, 1930, p. 361.
    * https://wiki.seg.org/wiki/Rock_physics#Voigt-Reuss-Hill_average
    """
    saturations = sum_saturations_to_one(saturations)
    return fluid_material(
        bulk_modulus=mix_bulk_moduli(fluids, saturations),
        density=mix_densities(fluids, saturations),
    )


def brie_fluid_mixing(fluids, fluid_saturations, gases, gas_saturations, exponent=3.0):
    """
    Fluid mixing law proposed by Brie et al. (1995), see also Borges et al. (2018).
    :param fluids: A list of liquids to mix with gas.
    :param fluid_saturations: List of saturation for each liquid.
    :param gases: The gases to mix into a liquid.
    :param gas_saturations: List of saturations for each gas. Total saturation
        (fluids + gases) should sum to 1.0.
    :param exponent: An empirical constant, defaults to 3.0.
    :return: Material representing the fluid mixed with gas.

    * Brie, A., et al. "Shear sonic interpretation in gas-bearing sands." SPE
        Annual Technical Conference and Exhibition. Society of Petroleum Engineers,
        1995.
    * Borges, F., and M. Landro. "Time-Lapse Separation Of Fluid And Pressure
        Effects With An Arbitrary Fluid Mixing Law." Fifth CO2 Geological Storage
        Workshop. Vol. 2018. No. 1. European Association of Geoscientists &
        Engineers, 2018.
    """
    if len(fluids) != len(fluid_saturations):
        raise ValueError("Mismatched lengths of fluid and fluid saturation input")
    if len(gases) != len(gas_saturations):
        raise ValueError("Mismatched lengths of gas and gas saturation input")

    if len(fluids) != 0:
        liquid = wood_fluid_mixing(fluids, fluid_saturations)
    if len(gases) != 0:
        gas = wood_fluid_mixing(gases, gas_saturations)

    if len(gases) == 0 and len(fluids) == 0:
        raise ValueError("Mixing is undefined for zero constituents")
    elif len(gases) == 0:
        return liquid
    elif len(fluids) == 0:
        return gas
    else:
        total_fluid = sum(fluid_saturations)
        bulk_mix = (liquid.bulk_modulus - gas.bulk_modulus) * (
            total_fluid
        ) ** exponent + gas.bulk_modulus
        density_mix = mix_densities([liquid, gas], [total_fluid, 1 - total_fluid])
        return fluid_material(density=density_mix, bulk_modulus=bulk_mix)


def fluid_substitution(dry_rock, mineral, fluid, porosity):
    """
    Performs fluid substitution according to gassmann theory
    see, for instance, "A Tutorial on Gassmann Fluid Substitution: Formulation,
    Algorithm and Matlab Code", Dhananjay Kumar, 2006,
    https://www.spgindia.org/geohorizon/jan_2006/dhananjay_paper.pdf

    :param dry_rock: A porous material without any pore-filling material.
    :param mineral: The material of the dry rock matrix.
    :param fluid: A fluid to saturate the dry rock with.
    :param porosity: The porosity of the dry rock as a ratio of vacuum.
    :returns: The resulting saturated material from injecting the dry rock with
        the given fluid.
    """

    @float_vectorize
    def gassman_equation(
        dry_bulk_modulus, mineral_bulk_modulus, fluid_bulk_modulus, porosity
    ):
        if porosity == 0:
            # If no porosity, no fluid injected, saturated material
            # is all mineral
            return mineral_bulk_modulus
        if dry_bulk_modulus == mineral_bulk_modulus:
            # Mineral is dry rock, so porosity can only be zero
            return mineral_bulk_modulus

        # To calculate saturated bulk modulus, we use the formulation from
        # Avseth et al.,"Quantitative Seismic Interpretation", page 17:
        # saturated_bulk_modulus / (mineral_bulk_modulus - saturated_bulk_modulus)
        # = dry_bulk_modulus / (mineral_bulk_modulus - dry_bulk_modulus) + fluid_bulk_modulus
        # / ((porosity * (mineral_bulk_modulus - fluid_bulk_modulus))
        difference_factor = dry_bulk_modulus / (
            mineral_bulk_modulus - dry_bulk_modulus
        ) + fluid_bulk_modulus / (
            (mineral_bulk_modulus - fluid_bulk_modulus) * porosity
        )
        # We now have saturated_bulk_modulus = difference_factor *
        # (mineral_bulk_modulus - saturated_bulk_modulus)

        if difference_factor < 0:
            raise ValueError(
                "moduli of materials are unsuitable for fluid substitution:\n"
                + f"fluid: {fluid_bulk_modulus}, mineral: {mineral_bulk_modulus},"
                + f"dry rock: {dry_bulk_modulus}"
            )

        return difference_factor / (1 + difference_factor) * mineral_bulk_modulus

    return Material(
        bulk_modulus=gassman_equation(
            dry_rock.bulk_modulus, mineral.bulk_modulus, fluid.bulk_modulus, porosity
        ),
        shear_modulus=dry_rock.shear_modulus,
        density=dry_rock.density + (porosity * fluid.density),
    )
