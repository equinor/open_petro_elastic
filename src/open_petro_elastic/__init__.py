from pkg_resources import DistributionNotFound, get_distribution

from .float_vectorize import float_vectorize

try:
    __version__ = get_distribution(__name__).version
except DistributionNotFound:
    __version__ = "0.0.0"

__all__ = [
    "float_vectorize",
]
