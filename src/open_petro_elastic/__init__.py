import importlib.metadata

from .float_vectorize import float_vectorize

__version__ = importlib.metadata.version(__package__ or __name__)

__all__ = [
    "float_vectorize",
]
