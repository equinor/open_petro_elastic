"""
Span-Wagner Coefficients

Coefficients defined in tables of Span & Wagner. Coefficients are divided into groups according to horizontal lines in
the table. This is also convenient when evaluating the Helmholtz energy expressions (phi functions).
"""

from numpy import array

""" Ideal gas part (table 27) """
a0 = array(
    [
        [
            8.37304456,
            -3.70454304,
            2.50000000,
            1.99427042,
            0.62105248,
            0.41195293,
            1.04028922,
            0.08327678,
        ]
    ]
)
theta0 = array(
    [
        [
            3.15163,
            6.11190,
            6.77708,
            11.32384,
            27.08792,
        ]
    ]
)

""" Residual part (table 31) """

# First group
n1 = array(
    [
        [
            0.38856823203161 * 1e0,
            0.29385475942740 * 1e1,
            -0.55867188534934 * 1e1,
            -0.76753199592477 * 1e0,
            0.31729005580416 * 1e0,
            0.54803315897767 * 1e0,
            0.12279411220335 * 1e0,
        ]
    ]
)
d1 = array([[1, 1, 1, 1, 2, 2, 3]])
t1 = array([[0.00, 0.75, 1.00, 2.00, 0.75, 2.00, 0.75]])

# Second group
n2 = array(
    [
        [
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
        ]
    ]
)
d2 = array(
    [[1, 2, 4, 5, 5, 5, 6, 6, 6, 1, 1, 4, 4, 4, 7, 8, 2, 3, 3, 5, 5, 6, 7, 8, 10, 4, 8]]
)
t2 = array(
    [
        [
            1.50,
            1.50,
            2.50,
            0.00,
            1.50,
            2.00,
            0.00,
            1.00,
            2.00,
            3.00,
            6.00,
            3.00,
            6.00,
            8.00,
            6.00,
            0.00,
            7.00,
            12.00,
            16.00,
            22.00,
            24.00,
            16.00,
            24.00,
            8.00,
            2.00,
            28.00,
            14.00,
        ]
    ]
)
c2 = array(
    [[1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 4, 4, 4, 4, 4, 4, 5, 6]]
)

# Third group
n3 = array(
    [
        [
            -0.21365488688320 * 1e3,
            0.26641569149272 * 1e5,
            -0.24027212204557 * 1e5,
            -0.28341603423999 * 1e3,
            0.21247284400179 * 1e3,
        ]
    ]
)
d3 = array([[2, 2, 2, 3, 3]])
t3 = array([[1.00, 0.00, 1.00, 3.00, 3.00]])
alpha3 = array([[25, 25, 25, 15, 20]])
beta3 = array([[325, 300, 300, 275, 275]])
gamma3 = array([[1.16, 1.19, 1.19, 1.25, 1.22]])
epsilon3 = array([[1.00, 1.00, 1.00, 1.00, 1.00]])

# Fourth group
n4 = array(
    [
        [
            -0.66642276540751 * 1e0,
            0.72608632349897 * 1e0,
            0.55068668612842 * 1e-1,
        ]
    ]
)
a4 = array([[3.500, 3.500, 3.000]])
b4 = array([[0.875, 0.925, 0.875]])
beta4 = array([[0.300, 0.300, 0.300]])
A4 = array([[0.700, 0.700, 0.700]])
B4 = array([[0.3, 0.3, 1.0]])
C4 = array([[10.0, 10.0, 12.5]])
D4 = array([[275, 275, 275]])
