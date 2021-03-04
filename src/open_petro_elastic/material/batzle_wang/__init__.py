"""
Calculation of pore fluid properties (oil, gas and brine) from

[1] Batzle, Michael, and Zhijing Wang. "Seismic properties of pore fluids."
Geophysics 57.11 (1992): 1396-1408.


see also https://petrowiki.spe.org/Pore_fluid_properties

"""

from .brine_properties import brine
from .hydro_carbon_gas import gas
from .oil_properties import dead_oil, live_oil, oil
from .water_properties import water

__all__ = ["dead_oil", "live_oil", "oil", "gas", "brine", "water"]
