from generators import materials, ratios
from hypothesis import assume, given
from predicates import between

from open_petro_elastic.material import hashin_shtrikman_walpole
from open_petro_elastic.material.sandstone import constant_cement, contact_cement


@given(
    ratios(min_value=0.2, max_value=0.49),
    ratios(min_value=0.1, max_value=0.3),
    ratios(),
    materials(),
    materials(),
)
def test_constant_cement(
    sand_porosity, contact_cement_porosity, cement_sand_ratio, sand, cement
):
    assume(sand_porosity > contact_cement_porosity)
    assume(cement.poisson_ratio >= 0)
    assume(sand.poisson_ratio >= 0)

    porosity = cement_sand_ratio * contact_cement_porosity
    sandstone = constant_cement(
        sand, cement, porosity, contact_cement_porosity, sand_porosity
    )

    dense_packing = hashin_shtrikman_walpole(
        cement, sand, sand_porosity - contact_cement_porosity
    )

    cc = contact_cement(cement, sand, contact_cement_porosity, sand_porosity)

    assert between(dense_packing.density, cc.density, sandstone.density)
    assert between(
        dense_packing.shear_modulus, cc.shear_modulus, sandstone.shear_modulus
    )
    assert between(dense_packing.bulk_modulus, cc.bulk_modulus, sandstone.bulk_modulus)
