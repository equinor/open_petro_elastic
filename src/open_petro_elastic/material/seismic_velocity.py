"""
Calculations derived from the usual equation for estimating primary
and secondary wave velocities in isotropic media:

Vp = sqrt((K+(4/3)u) / rho)

and

Vs = sqrt(u/rho)

where

  * Vp is the velocity of the primary (compressional) wave
  * Vs is the velocity of the secondary (shear) wave
  * K is the bulk moduli
  * u is the shear moduli


see also:

* "A Tutorial on Gassmann Fluid Substitution: Formulation, Algorithm and Matlab
Code", Dhananjay Kumar, 2006,
https://www.spgindia.org/geohorizon/jan_2006/dhananjay_paper.pdf

* https://petrowiki.org/Compressional_and_shear_velocities

"""


def calculate_shear_modulus(secondary_velocity, density):
    """
    Calculates shear modulus of a material from its density and velocity of
    secondary wave running through it.

    :param secondary_velocity: The velocity of secondary waves
        through the material in m/s.
    :param density: The density of the isotropic material in kg/m^3.
    :return: The shear modulus of the material
    """
    return secondary_velocity ** 2 * density


def calculate_bulk_modulus(primary_velocity, shear_modulus, density):
    """
    Calculates bulk modulus of a material from its density, shear modulus, and
    velocity of a primary wave running through it.
    :param primary_velocity: The velocity of primary waves through the material
    in m/s.
    :param density: The density of the isotropic material in kg/m^3.
    :param shear_modulus: The shear modulus of the material.
    :return: The bulk modulus of the material.
    """
    return primary_velocity ** 2 * density - 4 / 3 * shear_modulus


def calculate_primary_velocity(bulk_modulus, shear_modulus, density):
    """
    Calculates the velocity of primary waves in
    an isotropic material with the given moduli and density.
    :param shear_modulus: The shear modulus of the material.
    :param bulk_modulus: The bulk modulus of the material.
    :param density: The density of the isotropic material in kg/m^3.
    :return: The velocity of primary waves through the material
    in m/s.
    """
    return ((bulk_modulus + 4 / 3 * shear_modulus) / density) ** 0.5


def calculate_secondary_velocity(shear_modulus, density):
    """
    Calculates the velocity of primary waves in
    an isotropic material with the given shear modulus and density.
    :param shear_modulus: The shear modulus of the material.
    :param bulk_modulus: The bulk modulus of the material.
    :param density: The density of the isotropic material in kg/m^3.
    :return: The velocity of secondary waves through the material
    in m/s.
    """
    return (shear_modulus / density) ** 0.5
