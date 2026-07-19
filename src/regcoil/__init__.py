"""REGCOIL: A regularized current-potential method for computing stellarator coil shapes."""

from __future__ import annotations

from datetime import datetime
import logging
from typing import TextIO

__version__ = "0.1.0"

from . import _core
from .coil_surface import CoilSurface
from .fourier_surface import FourierSurface
from .log import log
from .plasma_surface import PlasmaSurface
from .regcoil import Regcoil, Solution
from .surface import Surface

__all__ = [
    "__version__",
    "_core",
    "log",
    "Surface",
    "FourierSurface",
    "PlasmaSurface",
    "CoilSurface",
    "Regcoil",
    "Solution",
]
