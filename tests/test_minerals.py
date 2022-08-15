import hypothesis.strategies as st
import pytest
from generators import constituents
from hypothesis import assume, given
from numpy.testing import assert_allclose
from predicates import assert_similar_material

from open_petro_elastic.config import Minerals
from open_petro_elastic.material import hashin_shtrikman_walpole


def test_no_constituents_raises_error():
    with pytest.raises(ValueError):
        Minerals([])


@given(constituents(), constituents(fraction=st.just(None)))
def test_mixed_is_hashin_shtrikman_walpole(sand, shale):
    assert_similar_material(
        Minerals([sand, shale]).as_mixture,
        hashin_shtrikman_walpole(sand.material, shale.material, sand.fraction),
    )


@given(constituents(fraction=st.just(None)))
def test_mixed_single_as_mixture_is_ident(shale):
    assert_similar_material(Minerals([shale]).as_mixture, shale.material)


@given(constituents(fraction=st.just(0.5)))
def test_minerals_non_one_sum(shale):
    with pytest.raises(ValueError, match="fraction"):
        Minerals([shale])


@given(constituents(), constituents(fraction=st.just(None)))
def test_mineral_from_config(sand, shale):
    sand_config = {
        "material": {
            "density": sand.material.density,
            "bulk_modulus": sand.material.bulk_modulus,
            "shear_modulus": sand.material.shear_modulus,
        },
        "fraction": sand.fraction,
    }
    shale_config = {
        "material": {
            "density": shale.material.density,
            "bulk_modulus": shale.material.bulk_modulus,
            "shear_modulus": shale.material.shear_modulus,
        },
        "fraction": shale.fraction,
    }

    assert_similar_material(
        Minerals([sand_config, shale_config]).as_mixture,
        Minerals([sand, shale]).as_mixture,
    )


@given(
    constituents(fraction=st.just(0.0)),
    constituents(fraction=st.just(0.0)),
    constituents(fraction=st.just(None)),
)
def test_mineral_zero_fraction_is_identity(sand, shale, quartz):
    assert_similar_material(
        Minerals([sand, shale, quartz]).as_mixture,
        quartz.material,
        atol=0.1,
    )


@given(
    constituents(),
    constituents(),
    constituents(fraction=st.just(None)),
)
def test_mineral_density_is_weighted_average(sand, shale, quartz):
    assume(sand.fraction + shale.fraction <= 1.0)
    minerals = Minerals([sand, shale, quartz])
    assert_allclose(
        minerals.as_mixture.density,
        sum(c.fraction * c.material.density for c in minerals),
    )
