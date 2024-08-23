import hypothesis.strategies as st
from hypothesis import assume

from open_petro_elastic.config import Constituent
from open_petro_elastic.material.material import Material


def positives(*args, **kwargs):
    defaults = {"allow_nan": False, "min_value": 0.1, "max_value": 10e8}
    defaults.update(kwargs)
    return st.floats(*args, **defaults)


def small_floats(*args, **kwargs):
    defaults = {"allow_nan": False, "width": 32, "allow_infinity": False}
    defaults.update(kwargs)
    return st.floats(*args, **defaults)


@st.composite
def materials(
    draw,
    bulk_modulus=positives(),
    shear_modulus=positives(),
    density=positives(),
):
    k = draw(bulk_modulus)
    mu = draw(shear_modulus)
    r = draw(density)
    try:
        m = Material(bulk_modulus=k, shear_modulus=mu, density=r)
    except ValueError:
        assume(False)
    return m


@st.composite
def fluids(draw, bulk_modulus=positives(), density=positives()):
    k = draw(bulk_modulus)
    r = draw(density)
    try:
        m = Material(bulk_modulus=k, shear_modulus=0, density=r)
    except ValueError:
        assume(False)
    return m


def ratios(*args, **kwargs):
    defaults = {"allow_nan": False, "min_value": 0.0, "max_value": 1.0}
    defaults.update(kwargs)
    return st.floats(*args, **defaults)


@st.composite
def constituents(draw, material=materials(), fraction=ratios()):
    return draw(st.builds(Constituent, material, fraction))
