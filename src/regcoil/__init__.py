"""REGCOIL: A regularized current-potential method for computing stellarator coil shapes."""

from __future__ import annotations

__version__ = "0.1.0"

from . import _core
from .coil_surface import CoilSurface
from .fourier_surface import FourierSurface
from .plasma_surface import PlasmaSurface
from .surface import Surface

__all__ = [
    "__version__",
    "_core",
    "Surface",
    "FourierSurface",
    "PlasmaSurface",
    "CoilSurface",
]
