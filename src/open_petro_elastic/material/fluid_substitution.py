from open_petro_elastic.float_vectorize import float_vectorize

from .material import Material


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
                + "dry rock: {dry_bulk_modulus}"
            )

        return difference_factor / (1 + difference_factor) * mineral_bulk_modulus

    return Material(
        bulk_modulus=gassman_equation(
            dry_rock.bulk_modulus, mineral.bulk_modulus, fluid.bulk_modulus, porosity
        ),
        shear_modulus=dry_rock.shear_modulus,
        density=dry_rock.density + (porosity * fluid.density),
    )
