"""
Package for calculating the residual helmholtz energy for CO2. The module uses sympy to evaluate the functions defined
in Span & Wagner [2], as well as its derivatives.
"""

import numpy as np
import sympy as sp

from . import coefficients


def residual_helmholtz_energy(delta_, tau_, diff_delta, diff_tau):
    """
    Equation 6.1 from Span & Wagner [2]. Calculates the residual helmholtz energy of co2 or its derivatives. tau_ and
    delta_ must have be numpy arrays of shape (N, 1). This allows for vectorization.

    :param delta_: Reduced density. Unit-less. numpy.ndarray with shape (N, 1)
    :param tau_: Inverse reduced temperature. Unit-less. numpy.ndarray with shape (N, 1)
    :param diff_delta: Degree of derivation wrt. delta. Integer.
    :param diff_tau: Degree of derivation wrt. tau. Integer.

    :return: Helmholtz free energy. Unit-less. numpy.ndarray with shape (N,)
    """
    _s1 = np.sum(
        _LAMBDIFIED_EXPRESSIONS[(s1, diff_delta, diff_tau)](tau_, delta_), axis=-1
    )
    _s2 = np.sum(
        _LAMBDIFIED_EXPRESSIONS[(s2, diff_delta, diff_tau)](tau_, delta_), axis=-1
    )
    _s3 = np.sum(
        _LAMBDIFIED_EXPRESSIONS[(s3, diff_delta, diff_tau)](tau_, delta_), axis=-1
    )
    _s4 = np.sum(
        _LAMBDIFIED_EXPRESSIONS[(s4, diff_delta, diff_tau)](tau_, delta_), axis=-1
    )

    return _s1 + _s2 + _s3 + _s4


# Define coefficients symbols. These should correspond one-to-one with variable names in the coefficient module
coeff_symbols = (
    n1,
    t1,
    d1,
    n2,
    d2,
    t2,
    c2,
    n3,
    d3,
    t3,
    alpha3,
    epsilon3,
    beta3,
    gamma3,
    n4,
    b4,
    a4,
    beta4,
    A4,
    B4,
    C4,
    D4,
) = sp.symbols(
    "n1 t1 d1 n2 d2 t2 c2 n3 d3 t3 alpha3 epsilon3 beta3 gamma3 n4 b4 a4 beta4 A4 B4 C4 D4"
)
tau, delta = sp.symbols("tau delta", real=True)
i = sp.symbols("i", integer=True)
s1 = n1 * delta**d1 * tau**t1
s2 = n2 * delta**d2 * tau**t2 * sp.exp(-(delta**c2))
s3 = (
    n3
    * delta**d3
    * tau**t3
    * sp.exp(-alpha3 * (delta - epsilon3) ** 2 - beta3 * (tau - gamma3) ** 2)
)

theta_expr = (1 - tau) + A4 * ((delta - 1) ** 2) ** (1 / (2 * beta4))
bigdelta_expr = theta_expr**2 + B4 * ((delta - 1) ** 2) ** a4
bigphi_expr = sp.exp(-C4 * (delta - 1) ** 2 - D4 * (tau - 1) ** 2)
s4 = n4 * bigdelta_expr**b4 * delta * bigphi_expr


def _lambdify(expr, diff_delta, diff_tau):
    diff = [delta] * diff_delta + [tau] * diff_tau
    if len(diff) > 0:
        expr = expr.diff(*diff)
    expr = expr.powsimp()
    coeff_vars = [getattr(coefficients, s.name) for s in coeff_symbols]
    _sympy_lambda = sp.utilities.lambdify(
        list(coeff_symbols) + [tau, delta],
        expr,
        modules=["numpy", {"DiracDelta": lambda x: x == 0}],
    )
    return lambda _tau, _delta: _sympy_lambda(*coeff_vars, _tau, _delta)


_LAMBDIFIED_EXPRESSIONS = {
    (e, dd, dt): _lambdify(e, dd, dt)
    for e in (s1, s2, s3, s4)
    for dd in (0, 1, 2)
    for dt in (0, 1, 2)
    if dd + dt <= 2
}
