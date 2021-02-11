from .fluid import fluid


def mix_densities(fluids, saturations):
    return sum(
        fluid.density * saturation for (fluid, saturation) in zip(fluids, saturations)
    )


def wood_fluid_mixing(fluids, saturations):
    """
    woods function for fluid mixing, also known as Reuss average.
    :param fluids: A list of liquids to mix with gas.
    :param saturations: List of saturation for each liquid (sums to 1.0).
    :return: Material representing the mixed fluids.

    * A. B. Wood, A Textbook of Sound, 1st ed. MacMillan, New York, 1930, p. 361.
    * https://wiki.seg.org/wiki/Rock_physics#Voigt-Reuss-Hill_average
    """
    bulk_mix = sum(saturations) / sum(
        saturation / fluid.bulk_modulus
        for (fluid, saturation) in zip(fluids, saturations)
    )
    return fluid(bulk_modulus=bulk_mix, density=mix_densities(fluids, saturations))


def brie_fluid_mixing(fluids, fluid_saturations, gas, exponent=3.0):
    """
    Fluid mixing law proposed by Brie et al. (1995), see also Borges et al. (2018).
    :param fluids: A list of liquids to mix with gas.
    :param saturations: List of saturation for each liquid ( sums to 1.0 minus saturation of gas).
    :param gas: The gas to mix into a liquid.
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
    liquid = wood_fluid_mixing(fluids, fluid_saturations)
    total_fluid = sum(fluid_saturations)
    bulk_mix = (liquid.bulk_modulus - gas.bulk_modulus) * (
        total_fluid
    ) ** exponent + gas.bulk_modulus
    density_mix = (
        mix_densities(fluids, fluid_saturations) + (1 - total_fluid) * gas.density
    )

    return fluid(density=density_mix, bulk_modulus=bulk_mix)
