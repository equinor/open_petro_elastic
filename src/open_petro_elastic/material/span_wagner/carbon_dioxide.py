"""
Equation 6.1 from Span & Wagner [2], and its derivatives have been computed using sympy and the following code:

First, define the kernels of the summations in the Helmholtz energy residual function of Table 32:
    >>> import sympy as sp
    ... (n1, t1, d1, n2, d2, t2, c2, n3, d3, t3, alpha3, epsilon3, beta3, gamma3, n4, b4, a4, beta4, A4, B4, C4, D4) =\
    ...     sp.symbols('n1 t1 d1 n2 d2 t2 c2 n3 d3 t3 alpha3 epsilon3 beta3 gamma3 n4 b4 a4 beta4 A4 B4 C4 D4')
    ... tau, delta = sp.symbols('tau delta', real=True)
    ... i = sp.symbols('i', integer=True)
    ... s1 = n1 * delta ** d1 * tau ** t1
    ... s2 = n2 * delta ** d2 * tau ** t2 * sp.exp(-delta ** c2)
    ... s3 = n3 * delta ** d3 * tau ** t3 * sp.exp(-alpha3 * (delta - epsilon3) **2 - beta3 * (tau - gamma3)**2)
    ...
    ... theta, bigdelta, bigphi = sp.symbols('theta bigdelta bigphi')
    ... theta_f, bigdelta_f, bigphi_f = sp.symbols('theta bigdelta bigphi', cls=sp.Function)
    ...
    ... theta_expr = (1 - tau) + A4 * ((delta - 1) ** 2) ** (1 / (2 * beta4))
    ... bigdelta_expr = theta_f(delta, tau) ** 2 + B4 * ((delta - 1) ** 2) ** a4
    ... bigphi_expr = sp.exp(-C4 * (delta - 1) ** 2 - D4 * (tau - 1) ** 2)
    ... s4 = n4 * bigdelta_f(delta, tau) ** b4 * delta * bigphi_f(delta, tau)
    ...
    ... # Print the expressions and its derivatives in a manner consistent with the code below
    ... # Define derivatives as symbols to enable printing to code
    ... dt_theta, dd_theta = sp.symbols('dt_theta dd_theta')
    ... dt_bigdelta, dd_bigdelta = sp.symbols('dt_bigdelta dd_bigdelta')
    ... dt_bigphi, dd_bigphi = sp.symbols('dt_bigphi dd_bigphi')
    ...
    ... dtt_theta, ddt_theta = sp.symbols('dtt_theta ddt_theta')
    ... dtt_bigdelta, ddt_bigdelta = sp.symbols('dtt_bigdelta ddt_bigdelta')
    ... dtt_bigphi, ddt_bigphi = sp.symbols('dtt_bigphi ddt_bigphi')
    ...
    ... ddd_theta = sp.symbols('ddd_theta')
    ... ddd_bigdelta = sp.symbols('ddd_bigdelta')
    ... ddd_bigphi = sp.symbols('ddd_bigphi')
    ...
    ... def print_expr(name, e, summation):
    ...     e = e.xreplace({
    ...         theta_f(delta, tau): theta,
    ...         bigdelta_f(delta, tau): bigdelta,
    ...         bigphi_f(delta, tau): bigphi,
    ...         theta_f(delta, tau).diff(tau): dt_theta,
    ...         theta_f(delta, tau).diff(delta): dd_theta,
    ...         bigdelta_f(delta, tau).diff(tau): dt_bigdelta,
    ...         bigdelta_f(delta, tau).diff(delta): dd_bigdelta,
    ...         bigphi_f(delta, tau).diff(tau): dt_bigphi,
    ...         bigphi_f(delta, tau).diff(delta): dd_bigphi,
    ...         theta_f(delta, tau).diff(tau, tau): dtt_theta,
    ...         theta_f(delta, tau).diff(delta, tau): ddt_theta,
    ...         bigdelta_f(delta, tau).diff(tau, tau): dtt_bigdelta,
    ...         bigdelta_f(delta, tau).diff(delta, tau): ddt_bigdelta,
    ...         bigphi_f(delta, tau).diff(tau, tau): dtt_bigphi,
    ...         bigphi_f(delta, tau).diff(delta, tau): ddt_bigphi,
    ...         theta_f(delta, tau).diff(delta, delta): ddd_theta,
    ...         bigdelta_f(delta, tau).diff(delta, delta): ddd_bigdelta,
    ...         bigphi_f(delta, tau).diff(delta, delta): ddd_bigphi,
    ...     })
    ...     e = e.powsimp()
    ...     kernel = str(e).replace("Abs", "abs")
    ...     if summation:
    ...         print(f'{name} = sum({kernel})')
    ...     else:
    ...         print(f'{name} = {kernel}')
    ...
    ... sep = '=' * 30
    ...
    ...
    ... # Base functions:
    ... print_expr('theta', theta_expr, False)
    ... print_expr('bigdelta', bigdelta_expr, False)
    ... print_expr('bigphi', bigphi_expr, False)
    ... for i, s in enumerate([s1, s2, s3, s4]):
    ...     print_expr(f's{i+1}', s, True)
    ... print(sep)
    ...
    ... # First derivatives
    ... print_expr('dt_theta', theta_expr.diff(tau), False)
    ... print_expr('dt_bigdelta', bigdelta_expr.diff(tau), False)
    ... print_expr('dt_bigphi', bigphi_expr.diff(tau), False)
    ... for i, s in enumerate([s1, s2, s3, s4]):
    ...     print_expr(f'dt_s{i+1}', s.diff(tau), True)
    ... print(sep)
    ... print_expr('dd_theta', theta_expr.diff(delta), False)
    ... print_expr('dd_bigdelta', bigdelta_expr.diff(delta), False)
    ... print_expr('dd_bigphi', bigphi_expr.diff(delta), False)
    ... for i, s in enumerate([s1, s2, s3, s4]):
    ...     print_expr(f'dd_s{i+1}', s.diff(delta), True)
    ... print(sep)
    ...
    ... # Second derivatives
    ... print_expr('dtt_theta', theta_expr.diff(tau, tau), False)
    ... print_expr('dtt_bigdelta', bigdelta_expr.diff(tau, tau), False)
    ... print_expr('dtt_bigphi', bigphi_expr.diff(tau, tau), False)
    ... for i, s in enumerate([s1, s2, s3, s4]):
    ...     print_expr(f'dtt_s{i+1}', s.diff(tau, tau), True)
    ... print(sep)
    ... print_expr('ddd_theta', theta_expr.diff(delta, delta), False)
    ... print_expr('ddd_bigdelta', bigdelta_expr.diff(delta, delta), False)
    ... print_expr('ddd_bigphi', bigphi_expr.diff(delta, delta), False)
    ... for i, s in enumerate([s1, s2, s3, s4]):
    ...     print_expr(f'ddd_s{i+1}', s.diff(delta, delta), True)
    ... print(sep)
    ... print_expr('ddt_theta', theta_expr.diff(delta, tau), False)
    ... print_expr('ddt_bigdelta', bigdelta_expr.diff(delta, tau), False)
    ... print_expr('ddt_bigphi', bigphi_expr.diff(delta, tau), False)
    ... for i, s in enumerate([s1, s2, s3, s4]):
    ...     print_expr(f'ddt_s{i+1}', s.diff(delta, tau), True)
    ... print(sep)

The expressions are directly equivalent to the expressions in the code below, with an exception being sympy's DiracDelta
function. This can, however, be directly replaced with a boolean expression.
"""
from numpy import array, log, exp, sum, abs, sign
from open_petro_elastic import float_vectorize

from open_petro_elastic.material.fluid import fluid_material as fluid
from open_petro_elastic.material.conversions import celsius_to_kelvin


"""
TODO:
    - Handle density (and other derived properties) close to phase transitions, critical point and triple point
    - Support vectorized variants of the below functions
"""


# Constants
R = 0.1889241 * 1e3  # J / kg K
T_critical = 304.1282  # K
rho_critical = 467.6  # kg / m^3
p_critical = 7.3773  # MPa
T_triple = 216.592  # K
p_triple = 0.51795  # MPa


"""
Span-Wagner Coefficients

Coefficients defined in Table 31 of Span & Wagner. Coefficients are divided into groups according to horizontal lines in
the table. This is also convenient when evaluating the Helmholtz energy expressions (phi functions).
"""


""" Ideal gas part """
a0 = array([
    8.37304456,
    -3.70454304,
    2.50000000,
    1.99427042,
    0.62105248,
    0.41195293,
    1.04028922,
    0.08327678,
])
theta0 = array([
    3.15163,
    6.11190,
    6.77708,
    11.32384,
    27.08792,
])


""" Residual part """
# First group
n1 = array([
    0.38856823203161 * 1e0,
    0.29385475942740 * 1e1,
    -0.55867188534934 * 1e1,
    -0.76753199592477 * 1e0,
    0.31729005580416 * 1e0,
    0.54803315897767 * 1e0,
    0.12279411220335 * 1e0,
])
d1 = array([1, 1, 1, 1, 2, 2, 3])
t1 = array([0.00, 0.75, 1.00, 2.00, 0.75, 2.00, 0.75])

# Second group
n2 = array([
    0.21658961543220 * 1e1,
    0.15841735109724 * 1e1,
    -0.23132705405503 * 1e0,
    0.58116916431436 * 1e-1,
    -0.55369137205382 * 1e0,
    0.48946615909422 * 1e0,
    -0.24275739843501 * 1e-1,
    0.62494790501678 * 1e-1,
    -0.12175860225246 * 1e0,
    -0.37055685270086 * 1e0,
    -0.16775879700426 * 1e-1,
    -0.11960736637987 * 1e0,
    -0.45619362508778 * 1e-1,
    0.35612789270346 * 1e-1,
    -0.74427727132052 * 1e-2,
    -0.17395704902432 * 1e-2,
    -0.21810121289527 * 1e-1,
    0.24332166559236 * 1e-1,
    -0.37440133423463 * 1e-1,
    0.143387157568781 * 1e0,
    -0.13491969083286 * 1e0,
    -0.23151225053480 * 1e-1,
    0.12363125492901 * 1e-1,
    0.21058321972940 * 1e-2,
    -0.33958519026368 * 1e-3,
    0.55993651771592 * 1e-2,
    -0.30335118055646 * 1e-3,
])
d2 = array([1, 2, 4, 5, 5, 5, 6, 6, 6, 1, 1, 4, 4, 4, 7, 8, 2, 3, 3, 5, 5, 6, 7, 8, 10, 4, 8])
t2 = array([1.50, 1.50, 2.50, 0.00, 1.50, 2.00, 0.00, 1.00, 2.00, 3.00, 6.00, 3.00, 6.00, 8.00, 6.00, 0.00, 7.00, 12.00,
            16.00, 22.00, 24.00, 16.00, 24.00, 8.00, 2.00, 28.00, 14.00])
c2 = array([1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 4, 4, 4, 4, 4, 4, 5, 6])

# Third group
n3 = array([
    -0.21365488688320 * 1e3,
    0.26641569149272 * 1e5,
    -0.24027212204557 * 1e5,
    -0.28341603423999 * 1e3,
    0.21247284400179 * 1e3,
])
d3 = array([2, 2, 2, 3, 3])
t3 = array([1.00, 0.00, 1.00, 3.00, 3.00])
alpha3 = array([25, 25, 25, 15, 20])
beta3 = array([325, 300, 300, 275, 275])
gamma3 = array([1.16, 1.19, 1.19, 1.25, 1.22])
epsilon3 = array([1.00, 1.00, 1.00, 1.00, 1.00])

# Fourth group
n4 = array([
    -0.66642276540751 * 1e0,
    0.72608632349897 * 1e0,
    0.55068668612842 * 1e-1,
])
a4 = array([3.500, 3.500, 3.000])
b4 = array([0.875, 0.925, 0.875])
beta4 = array([0.300, 0.300, 0.300])
A4 = array([0.700, 0.700, 0.700])
B4 = array([0.3, 0.3, 1.0])
C4 = array([10.0, 10.0, 12.5])
D4 = array([275, 275, 275])


def phi(delta, tau, dd, dt):
    """
    Helmholtz energy as defined by equation 6.1 in Span & Wagner [2]

    :param delta: Reduced density, unit-less. That is, density / rho_critical.
    :param tau: Inverse reduced temperature, unit-less. That is, T_critical / (absolute) temperature
    :param dd: Degree of derivation wrt. delta. Integer between 0 and 2.
    :param dt: Degree of derivation wrt. tau. Integer between 0 and 2, as long as (dt + dd < 3)
    """
    return phi_0(delta, tau, dd, dt) + phi_r(delta, tau, dd, dt)


def phi_0(delta, tau, dd, dt):
    """
    Helmholtz energy from ideal gas behavior as defined by equation 2.3 in Span & Wagner [2]. See function phi for
    argument documentation.
    """
    if dt == dd == 0:
        _sum = sum(a0[3:] * log(1 - exp(-theta0 * tau)))
        return log(delta) + a0[0] + a0[1] * tau + a0[2] * log(tau) + _sum
    elif dt == 1 and dd == 0:
        _sum = sum(a0[3:] * theta0 * ((1 - exp(-theta0 * tau)) ** -1 - 1))
        return a0[1] + a0[2] / tau + _sum
    elif dt == 0 and dd == 1:
        return 1 / delta
    elif dt == 2 and dd == 0:
        _sum = sum(a0[3:] * theta0 ** 2 * exp(-theta0 * tau) * (1 - exp(-theta0 * tau)) ** -2)
        return -a0[2] / tau ** 2 - _sum
    elif dt == 0 and dd == 2:
        return - 1 / delta ** 2
    elif dt == 1 and dd == 1:
        return 0
    else:
        raise NotImplementedError


def phi_r(delta, tau, dd, dt):
    """
    Residual part of Helmholtz energy as defined by the equation in Table 32 of Span & Wagner [2]. See function phi for
    argument documentation.
    """
    # tau == 1.0 or delta == 1.0 leads to numerically invalid results. The values are nudged to avoid nan output.
    if tau == 1.0:
        tau -= 1e-15
    if delta == 1.0:
        delta -= 1e-15

    theta = A4*abs(delta - 1)**(1/beta4) - tau + 1
    bigdelta = B4*abs(delta - 1)**(2*a4) + theta**2
    bigphi = exp(-C4*(delta - 1)**2 - D4*(tau - 1)**2)
    if dt == dd == 0:
        s1 = sum(delta**d1*n1*tau**t1)
        s2 = sum(delta**d2*n2*tau**t2*exp(-delta**c2))
        s3 = sum(delta**d3*n3*tau**t3*exp(-alpha3*(delta - epsilon3)**2 - beta3*(-gamma3 + tau)**2))
        s4 = sum(bigdelta**b4*bigphi*delta*n4)
        return s1 + s2 + s3 + s4
    elif dt == 1 and dd == 0:
        dt_theta = -1
        dt_bigdelta = 2*dt_theta*theta
        dt_bigphi = -D4*(2*tau - 2)*exp(-C4*(delta - 1)**2 - D4*(tau - 1)**2)
        dt_s1 = sum(delta**d1*n1*t1*tau**(t1 - 1))
        dt_s2 = sum(delta**d2*n2*t2*tau**(t2 - 1)*exp(-delta**c2))
        dt_s3 = sum(-beta3*delta**d3*n3*tau**t3*(-2*gamma3 + 2*tau)*exp(-alpha3*(delta - epsilon3)**2 - beta3*(
                -gamma3 + tau)**2) + delta**d3*n3*t3*tau**(t3 - 1)*exp(-alpha3*(delta - epsilon3)**2 - beta3*(-gamma3 +
                                                                                                              tau)**2))
        dt_s4 = sum(b4*bigdelta**(b4 - 1)*bigphi*delta*dt_bigdelta*n4 + bigdelta**b4*delta*dt_bigphi*n4)
        return dt_s1 + dt_s2 + dt_s3 + dt_s4
    elif dt == 0 and dd == 1:
        dd_theta = A4 * abs(delta - 1) ** (-1 + 1 / beta4) * sign(delta - 1) / beta4
        dd_bigdelta = 2 * B4 * a4 * abs(delta - 1) ** (2 * a4 - 1) * sign(delta - 1) + 2 * dd_theta * theta
        dd_bigphi = -C4 * (2 * delta - 2) * exp(-C4 * (delta - 1) ** 2 - D4 * (tau - 1) ** 2)
        dd_s1 = sum(d1 * delta ** (d1 - 1) * n1 * tau ** t1)
        dd_s2 = sum(-c2 * delta ** (c2 + d2 - 1) * n2 * tau ** t2 * exp(-delta ** c2) + d2 * delta ** (
                    d2 - 1) * n2 * tau ** t2 * exp(-delta ** c2))
        dd_s3 = sum(-alpha3 * delta ** d3 * n3 * tau ** t3 * (2 * delta - 2 * epsilon3) * exp(
            -alpha3 * (delta - epsilon3) ** 2 - beta3 * (-gamma3 + tau) ** 2) + d3 * delta ** (
                                d3 - 1) * n3 * tau ** t3 * exp(
            -alpha3 * (delta - epsilon3) ** 2 - beta3 * (-gamma3 + tau) ** 2))
        dd_s4 = sum(b4 * bigdelta ** (
                    b4 - 1) * bigphi * dd_bigdelta * delta * n4 + bigdelta ** b4 * bigphi * n4 + bigdelta ** b4 *
                    dd_bigphi * delta * n4)
        return dd_s1 + dd_s2 + dd_s3 + dd_s4
    else:
        dt_theta = -1
        dt_bigdelta = 2*dt_theta*theta
        dt_bigphi = -D4*(2*tau - 2)*exp(-C4*(delta - 1)**2 - D4*(tau - 1)**2)
        dd_theta = A4 * abs(delta - 1) ** (-1 + 1 / beta4) * sign(delta - 1) / beta4
        dd_bigdelta = 2 * B4 * a4 * abs(delta - 1) ** (2 * a4 - 1) * sign(delta - 1) + 2 * dd_theta * theta
        dd_bigphi = -C4 * (2 * delta - 2) * exp(-C4 * (delta - 1) ** 2 - D4 * (tau - 1) ** 2)
        if dt == 2 and dd == 0:
            dtt_theta = 0
            dtt_bigdelta = 2 * dt_theta ** 2 + 2 * dtt_theta * theta
            dtt_bigphi = 2 * D4 * (2 * D4 * (tau - 1) ** 2 - 1) * exp(-C4 * (delta - 1) ** 2 - D4 * (tau - 1) ** 2)
            dtt_s1 = sum(delta ** d1 * n1 * t1 * tau ** (t1 - 2) * (t1 - 1))
            dtt_s2 = sum(delta ** d2 * n2 * t2 * tau ** (t2 - 2) * (t2 - 1) * exp(-delta ** c2))
            dtt_s3 = sum(delta ** d3 * n3 * tau ** t3 * (4 * beta3 * t3 * (gamma3 - tau) / tau + 2 * beta3 * (
                        2 * beta3 * (gamma3 - tau) ** 2 - 1) + t3 * (t3 - 1) / tau ** 2) * exp(
                -alpha3 * (delta - epsilon3) ** 2 - beta3 * (-gamma3 + tau) ** 2))
            dtt_s4 = sum(bigdelta ** b4 * delta * n4 * (b4 * bigphi * (
                        b4 * dt_bigdelta ** 2 / bigdelta + dtt_bigdelta - dt_bigdelta ** 2 / bigdelta) / bigdelta + 2 *
                                                        b4 * dt_bigdelta * dt_bigphi / bigdelta + dtt_bigphi))
            return dtt_s1 + dtt_s2 + dtt_s3 + dtt_s4
        elif dt == 0 and dd == 2:
            ddd_theta = A4 * (
                        2 * (delta == 1) / abs(delta - 1) - sign(delta - 1) ** 2 / (delta - 1) ** 2 + sign(
                    delta - 1) ** 2 / (beta4 * (delta - 1) ** 2)) * abs(delta - 1) ** (1 / beta4) / beta4
            ddd_bigdelta = 4 * B4 * a4 ** 2 * abs(delta - 1) ** (2 * a4) * sign(delta - 1) ** 2 / (
                        delta - 1) ** 2 + 4 * B4 * a4 * abs(delta - 1) ** (2 * a4 - 1) * (delta == 1) - 2 * B4 * a4 *\
                           abs(delta - 1) ** (2 * a4) * sign(delta - 1) ** 2 / (
                                       delta - 1) ** 2 + 2 * dd_theta ** 2 + 2 * ddd_theta * theta
            ddd_bigphi = 2 * C4 * (2 * C4 * (delta - 1) ** 2 - 1) * exp(-C4 * (delta - 1) ** 2 - D4 * (tau - 1) ** 2)
            ddd_s1 = sum(d1 * delta ** (d1 - 2) * n1 * tau ** t1 * (d1 - 1))
            ddd_s2 = sum(delta ** (d2 - 2) * n2 * tau ** t2 * (
                        -2 * c2 * d2 * delta ** c2 + c2 * delta ** c2 * (c2 * delta ** c2 - c2 + 1) + d2 * (
                            d2 - 1)) * exp(-delta ** c2))
            ddd_s3 = sum(delta ** d3 * n3 * tau ** t3 * (-4 * alpha3 * d3 * (delta - epsilon3) / delta + 2 * alpha3 * (
                        2 * alpha3 * (delta - epsilon3) ** 2 - 1) + d3 * (d3 - 1) / delta ** 2) * exp(
                -alpha3 * (delta - epsilon3) ** 2 - beta3 * (-gamma3 + tau) ** 2))
            ddd_s4 = sum(bigdelta ** b4 * n4 * (2 * b4 * bigphi * dd_bigdelta / bigdelta + b4 * bigphi * delta * (
                        b4 * dd_bigdelta ** 2 / bigdelta + ddd_bigdelta - dd_bigdelta ** 2 / bigdelta) / bigdelta + 2 *
                                                b4 * dd_bigdelta * dd_bigphi * delta / bigdelta + 2 * dd_bigphi +
                                                ddd_bigphi * delta))
            return ddd_s1 + ddd_s2 + ddd_s3 + ddd_s4
        elif dt == 1 and dd == 1:
            ddt_theta = 0
            ddt_bigdelta = 2 * dd_theta * dt_theta + 2 * ddt_theta * theta
            ddt_bigphi = 4 * C4 * D4 * (delta - 1) * (tau - 1) * exp(-C4 * (delta - 1) ** 2 - D4 * (tau - 1) ** 2)
            ddt_s1 = sum(d1 * delta ** (d1 - 1) * n1 * t1 * tau ** (t1 - 1))
            ddt_s2 = sum(delta ** (d2 - 1) * n2 * t2 * tau ** (t2 - 1) * (-c2 * delta ** c2 + d2) * exp(-delta ** c2))
            ddt_s3 = sum(delta ** d3 * n3 * tau ** t3 * (
                        -4 * alpha3 * beta3 * (delta - epsilon3) * (gamma3 - tau) - 2 * alpha3 * t3 * (
                            delta - epsilon3) / tau + 2 * beta3 * d3 * (gamma3 - tau) / delta + d3 * t3 / (
                                    delta * tau)) * exp(
                -alpha3 * (delta - epsilon3) ** 2 - beta3 * (-gamma3 + tau) ** 2))
            ddt_s4 = sum(bigdelta ** b4 * n4 * (
                        b4 ** 2 * bigphi * dd_bigdelta * delta * dt_bigdelta / bigdelta ** 2 + b4 * bigphi *
                        ddt_bigdelta * delta / bigdelta + b4 * bigphi * dt_bigdelta / bigdelta + b4 * dd_bigdelta *
                        delta * dt_bigphi / bigdelta + b4 * dd_bigphi * delta * dt_bigdelta / bigdelta - b4 * bigphi *
                        dd_bigdelta * delta * dt_bigdelta / bigdelta ** 2 + ddt_bigphi * delta + dt_bigphi))
            return ddt_s1 + ddt_s2 + ddt_s3 + ddt_s4
        else:
            raise NotImplementedError


@float_vectorize
def carbon_dioxide_pressure(absolute_temperature, density, d_density=0, d_temperature=0):
    """
    CO2 pressure (MPa) as given by Table 3 of Span & Wagner [2]

    :param absolute_temperature: Temperature in K.
    :param density: CO2 density (kg / m^3).
    :param d_density: Degree of derivation wrt. density.
    :param d_temperature: Degree of derivation wrt. temperature.
    """
    # Table 3, SW
    tau = T_critical / absolute_temperature
    delta = density / rho_critical
    if d_temperature != 0:
        raise NotImplementedError
    if d_density == 0:
        return density * R * absolute_temperature * (1 + delta * phi_r(delta, tau, 1, 0)) / 1e6
    elif d_density == 1:
        first = 2 * delta * phi_r(delta, tau, 1, 0)
        second = density * delta * phi_r(delta, tau, 2, 0)
        return absolute_temperature * R * (1 + first + second) / 1e6


def saturated_liquid_density(absolute_temperature):
    """
    Saturated liquid density as defined by equation 3.14 of Span & Wagner [2]

    :param absolute_temperature: Absolute temperature in K. Should satisfy T_triple < absolute_temperature < T_critical
    """
    _a1 = 1.9245108
    _a2 = -0.62385555
    _a3 = -0.32731127
    _a4 = 0.39245142
    _t = 1 - absolute_temperature / T_critical
    inner = _a1 * _t ** 0.34 + _a2 * _t ** 0.5 + _a3 * _t ** (10 / 6) + _a4 * _t ** (11 / 6)
    return rho_critical * exp(inner)


def saturated_vapor_density(absolute_temperature):
    """
    Saturated vapor density as defined by equation 3.15 of Span & Wagner

    :param absolute_temperature: Absolute temperature in K. Should satisfy T_triple < absolute_temperature < T_critical
    """
    # Assert temp < critical
    _a1 = -1.7074879
    _a2 = -0.82274670
    _a3 = -4.6008549
    _a4 = -10.111178
    _a5 = -29.742252
    _t = 1 - absolute_temperature / T_critical
    inner = _a1 * _t ** 0.34 + _a2 * _t ** 0.5 + _a3 * _t + _a4 * _t ** (7 / 3) + _a5 * _t ** (14 / 3)
    return rho_critical * exp(inner)


def sublimation_pressure(absolute_temperature):
    """
    Sublimation pressure as defined by equation 3.12 of Span & Wagner [2]

    :param absolute_temperature: Absolute temperature in K. Should satisfy absolute_temperature < T_triple
    """
    _a1 = -14.740846
    _a2 = 2.4327015
    _a3 = -5.3061778
    _t = 1 - absolute_temperature / T_triple
    inner = _a1 * _t + _a2 * _t ** 1.9 + _a3 * _t ** 2.9
    return p_triple * exp(inner / (1 - _t))


def vapor_pressure(absolute_temperature):
    """
    Vapor pressure as defined by equation 3.13 of Span & Wagner [2]

    :param absolute_temperature: Absolute temperature in K. Should satisfy T_triple < absolute_temperature < T_critical
    """
    _a1 = -7.0602087
    _a2 = 1.9391218
    _a3 = -1.6463597
    _a4 = -3.2995634
    _t = 1 - absolute_temperature / T_critical
    inner = _a1 * _t ** 1.0 + _a2 * _t ** 1.5 + _a3 * _t ** 2.0 + _a4 * _t ** 4.0
    return p_critical * exp(inner / (1 - _t))


def melting_pressure(absolute_temperature):
    """
    Melting pressure as defined by equation 3.10 of Span & Wagner [2]

    :param absolute_temperature: Absolute temperature in K. Should satisfy T_triple < absolute_temperature
    """
    _a1 = 1955.5390
    _a2 = 2055.4593
    _t = absolute_temperature / T_triple - 1
    return p_triple * (1 + _a1 * _t + _a2 * _t ** 2)


@float_vectorize
def carbon_dioxide_density(absolute_temperature, pressure, force_vapor='auto'):
    """
    Density of carbon dioxide. Found solving the Pressure equation of Table 3 in Span & Wagner [2] numerically for
    density. To ensure a single solution, the phase of the liquid must first be determined.

    :param absolute_temperature: Absolute temperature (K).
    :param pressure: Pressure (MPa).
    :param force_vapor: Boolean or 'auto'. If 'auto', the phase of the fluid is automatically determined. However, along
        the vaporization line (assuming T_triple < absolute_temperature < T_critical), the fluid is in two-phase
        equilibrium and the phase cannot be uniquely determined. If force_vapor is set to True, vapor phase is always
        selected, if False, liquid phase is selected. Outide the temperature bounds, this argument has no effect. This
        argument should only be used close to the vaporization boundary, otherwise the behavior might not be as
        expected.
    :return: Density (kg / m^3)
    """
    import scipy.optimize
    bounds = [0.1, 1500.0]  # TODO: these are set through inspection. Should check publication for bounds
    if absolute_temperature < T_triple:
        bounds[1] = saturated_vapor_density(absolute_temperature)
    elif absolute_temperature < T_critical:
        if (force_vapor == 'auto' and pressure < vapor_pressure(absolute_temperature)) or force_vapor is True:
            bounds[1] = saturated_vapor_density(absolute_temperature)
        else:
            bounds[0] = saturated_liquid_density(absolute_temperature)

    # Extend bounds slightly to ensure toms748 can converge
    bounds = [bounds[0] * 0.95, bounds[1] * 1.05]

    opt = scipy.optimize.root_scalar(
        lambda x: carbon_dioxide_pressure(absolute_temperature, x) - pressure,
        # fprime=lambda x: carbon_dioxide_pressure(absolute_temperature, x, 1, 0),
        method='toms748',
        bracket=bounds,
        x0=sum(bounds) / 2,
    )
    return opt.root


def carbon_dioxide_bulk_modulus(absolute_temperature, pressure):
    raise NotImplementedError


def carbon_dioxide(temperature, pressure):
    """
    :param temperature: Temperature in Celsius. Valid range: 216 K - 1100 K
    :param pressure: Pressure in MPa. Valid range: 0 - 800
    """
    t = celsius_to_kelvin(temperature)
    return fluid(
        density=carbon_dioxide_density(t, pressure),
        bulk_modulus=carbon_dioxide_bulk_modulus(t, pressure)
    )
