"""REGCOIL: A regularized current-potential method for computing stellarator coil shapes."""

from __future__ import annotations

__version__ = "0.1.0"

from . import _core
from .coil_surface import CoilSurface
from .fourier_surface import FourierSurface
from .log import log
from .plasma_surface import PlasmaSurface
from .regcoil import Regcoil, Solution
from .surface import Surface
from ._core import omp_max_threads

__all__ = [
    "__version__",
    "_core",
    "omp_max_threads",
    "log",
    "Surface",
    "FourierSurface",
    "PlasmaSurface",
    "CoilSurface",
    "Regcoil",
    "Solution",
]
