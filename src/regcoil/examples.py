"""Bundled example datasets (VMEC equilibria, BNORM, and winding-surface files).

The data files live in the source tree under ``equilibria/`` and are installed
into the package as ``regcoil/equilibria/`` (see ``src/regcoil/meson.build``).
Reach them at runtime through the :data:`examples` registry::

    import regcoil

    ds = regcoil.examples("W7-X")
    plasma = regcoil.PlasmaSurface.from_wout(ds.wout, ntheta=128, nzeta=128)
    plasma.set_bnormal_from_bnorm_file(ds.bnorm)
    coil = regcoil.CoilSurface.from_nescin(ds.nescin, nfp=plasma.nfp)

Lookups are case-insensitive. The accepted names (including aliases) are
returned by ``regcoil.examples.available()``.
"""

from __future__ import annotations

from dataclasses import dataclass
from importlib.resources import files
from pathlib import Path

__all__ = ["ExampleDataset", "examples"]


@dataclass(frozen=True)
class ExampleDataset:
    """Filesystem paths to the files that make up one example configuration.

    ``name`` is the canonical dataset id; the other fields are absolute
    :class:`~pathlib.Path` objects that can be passed straight to the
    ``from_wout`` / ``set_bnormal_from_bnorm_file`` / ``from_nescin`` loaders.
    """

    name: str
    wout: Path
    bnorm: Path
    nescin: Path


# Canonical dataset id -> (wout, bnorm, nescin) filenames within equilibria/.
_REGISTRY: dict[str, tuple[str, str, str]] = {
    "li383_1.4m": (
        "wout_li383_1.4m.nc",
        "bnorm.li383_1.4m",
        "nescin.li383_realWindingSurface",
    ),
    "d23p4_tm": (
        "wout_d23p4_tm.nc",
        "bnorm.d23p4_tm",
        "nescin.w7x_winding_surface_from_Drevlak",
    ),
}

# Lowercased lookup key -> canonical dataset id. Every canonical id is also a
# valid key so a dataset can always be requested by its own name.
_ALIASES: dict[str, str] = {
    "ncsx": "li383_1.4m",
    "li383": "li383_1.4m",
    "li383_1.4m": "li383_1.4m",
    "w7-x": "d23p4_tm",
    "w7x": "d23p4_tm",
    "d23p4_tm": "d23p4_tm",
}

# Human-facing names shown by available() and in error messages, in the order
# a user is most likely to reach for them.
_DISPLAY_NAMES = ["NCSX", "li383", "li383_1.4m", "W7-X", "d23p4_tm"]


def _resolve(filename: str) -> Path:
    """Return the installed filesystem path of one bundled data file.

    meson-python installs the package (and this data) as ordinary files on
    disk, never inside a zip, so the resource is always a real path.
    """
    resource = files("regcoil").joinpath("equilibria", filename)
    if not resource.is_file():
        raise FileNotFoundError(
            f"Bundled example data file {filename!r} was not found at "
            f"{resource}. If you are running from a source checkout, the "
            f"package may need to be (re)installed so meson copies the "
            f"equilibria/ data into it."
        )
    return Path(str(resource))


class _Examples:
    """Callable registry of bundled example datasets.

    Call it with a dataset name or alias to get an :class:`ExampleDataset`;
    call :meth:`available` to list the accepted names.
    """

    def __call__(self, name: str) -> ExampleDataset:
        key = _ALIASES.get(name.strip().lower())
        if key is None:
            raise KeyError(
                f"Unknown example dataset {name!r}. "
                f"Available: {', '.join(_DISPLAY_NAMES)}."
            )
        wout, bnorm, nescin = _REGISTRY[key]
        return ExampleDataset(
            name=key,
            wout=_resolve(wout),
            bnorm=_resolve(bnorm),
            nescin=_resolve(nescin),
        )

    @staticmethod
    def available() -> list[str]:
        """Human-facing dataset names and aliases accepted by the registry."""
        return list(_DISPLAY_NAMES)

    def __repr__(self) -> str:  # pragma: no cover - cosmetic
        return f"<regcoil examples: {', '.join(_DISPLAY_NAMES)}>"


examples = _Examples()
