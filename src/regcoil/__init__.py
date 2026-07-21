"""REGCOIL: A regularized current-potential method for computing stellarator coil shapes."""

from __future__ import annotations

__version__ = "0.1.0"

from . import _core
from . import cut
from . import plot
from .coil_surface import CoilSurface
from .cut import CutCoils
from .examples import ExampleDataset, examples
from .fourier_surface import FourierSurface
from .log import log
from .plasma_surface import PlasmaSurface
from .regcoil import Regcoil, Solution, SolutionScan
from .reparameterize import CurvatureWeighted, ThetaMap, UniformArclength
from .surface import Surface
from ._core import omp_max_threads
from ._serialize import load, save

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
    "SolutionScan",
    "CutCoils",
    "UniformArclength",
    "CurvatureWeighted",
    "ThetaMap",
    "ExampleDataset",
    "examples",
    "plot",
    "cut",
    "save",
    "load",
]
