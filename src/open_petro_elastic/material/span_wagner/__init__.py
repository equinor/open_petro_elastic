"""
Calculation of CO2 properties from

[2] Span, Roland, and Wolfgang Wagner. "A new equation of state for carbon dioxide covering the fluid region from the
triple‐point temperature to 1100 K at pressures up to 800 MPa." Journal of physical and chemical reference data 25.6
(1996): 1509-1596.
"""

from .carbon_dioxide import carbon_dioxide

__all__ = ["carbon_dioxide"]
