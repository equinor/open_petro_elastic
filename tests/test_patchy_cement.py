from generators import materials, positives, ratios
from hypothesis import assume, given

from open_petro_elastic.material.sandstone import patchy_cement


@given(
    ratios(min_value=0.2, max_value=0.49),
    ratios(min_value=0.1, max_value=0.3),
    ratios(min_value=0.1, max_value=0.3),
    positives(min_value=1.0),
    positives(min_value=1.0),
    ratios(max_value=0.9),
    materials(),
    materials(),
)
def test_patchy_cement(
    sand_porosity,
    cemented_sand_porosity1,
    cemented_sand_porosity2,
    pressure1,
    pressure2,
    cement_sand_ratio,
    sand,
    cement,
):
    lower_bound_cemented_sand_porosity = min(
        cemented_sand_porosity1, cemented_sand_porosity2
    )
    cemented_sand_porosity = max(cemented_sand_porosity1, cemented_sand_porosity2)
    lower_bound_pressure = min(pressure1, pressure2)
    pressure = max(pressure1, pressure2)
    assume(pressure > lower_bound_pressure)
    assume(sand_porosity > cemented_sand_porosity)
    porosity = cement_sand_ratio * lower_bound_cemented_sand_porosity
    patchy_cement(
        sand,
        cement,
        porosity,
        cemented_sand_porosity,
        lower_bound_cemented_sand_porosity,
        pressure,
        lower_bound_pressure,
        sand_porosity,
    )
