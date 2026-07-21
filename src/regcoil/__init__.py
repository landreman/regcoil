"""REGCOIL: A regularized current-potential method for computing stellarator coil shapes."""

from __future__ import annotations

# Single source of truth for the version is [project.version] in pyproject.toml;
# this reads it back out of the installed distribution metadata so the two can
# never drift.  regcoil always has to be installed (it carries the compiled
# _core extension), so there is no source-tree case to fall back on.
from importlib.metadata import version as _version

__version__ = _version("regcoil")

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
