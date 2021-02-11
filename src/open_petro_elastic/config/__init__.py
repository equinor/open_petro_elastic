"""
The open_petro_elastic.config module contains the internal representation of
the configuration file. The parsing of the configuration file goes through
three representations 1) YAML, 2) dictionary of yaml layout 3) Internal representation
as objects built from each section of the config file.

Each of the objects can be built from dictionaries, but it is possible to
construct directly. Generally, the constructor takes the unpacked
corresponding section of the dictionary, ie.  Fluids(**config["fluids"]) where
config is the dictionary from step 2.

"""
from .coefficients import Coefficients
from .constituent import Constituent
from .depth_trend import DepthTrend
from .dry_rock import DryRock
from .fluids import Fluids
from .minerals import Minerals
from .pressure import Pressure
from .pressure_dependency import PressureDependency

__all__ = [
    "DryRock",
    "Pressure",
    "Coefficients",
    "PressureDependency",
    "DepthTrend",
    "Minerals",
    "Constituent",
    "Fluids",
]
