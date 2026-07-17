"""regcoil – regularized current-potential method for stellarator coil shapes."""

from __future__ import annotations

__version__ = "0.1.0"

from . import _core
from ._core import RegcoilProblem

__all__ = ["__version__", "_core", "RegcoilProblem"]
