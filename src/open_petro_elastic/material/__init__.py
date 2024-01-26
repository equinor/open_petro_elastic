"""
The materials submodule contains the models for materials used
for open_petro_elastic and can be used as a standalone library.

Sentral to the library is the Material class which carries the
elastic properties of any given material

    >>> from open_petro_elastic.material.sandstone import hertz_mindlin
    >>> from open_petro_elastic.material import Material
    >>> mineral = Material(
    ...        density=2.75,
    ...        bulk_modulus=20000,
    ...        shear_modulus=7,
    ... )
    >>> hertz_mindlin(mineral, 0.4, 1.0)
    Material(shear_modulus=3.8..., bulk_modulus=3.179...,...density=1.65)


"""

from .fluid import brie_fluid_mixing, fluid_substitution, wood_fluid_mixing
from .hashin_shtrikman import hashin_shtrikman_average, hashin_shtrikman_walpole
from .material import Material, polyval_material

__all__ = [
    "Material",
    "hashin_shtrikman_average",
    "hashin_shtrikman_walpole",
    "fluid_substitution",
    "brie_fluid_mixing",
    "wood_fluid_mixing",
    "polyval_material",
]
