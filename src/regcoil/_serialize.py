"""Save/load: object serialization to NetCDF-4 via `h5netcdf` (ADR-028,
Phase 9). See docs/migration/API.md#saving-and-loading for the schema.

Files are grouped by object (`/plasma`, `/coil`, `/problem`, `/solutions`),
each tagged with a `_class` attribute; the root carries `format_version`.
Only the source of truth, plus derived quantities whose recompute needs the
Fortran kernel or the big operators, are stored -- everything else is
rebuilt on load from the surfaces and mode amplitudes.
"""

from __future__ import annotations

import numpy as np

FORMAT_VERSION = 1


def _to_bool_attr(value):
    return np.int8(1 if value else 0)


def _write_fourier_surface(group, surf, dim_prefix):
    mnmax_dim = f"mnmax_{dim_prefix}"
    group.dimensions[mnmax_dim] = surf.mnmax
    for name in ("xm", "xn", "rmnc", "rmns", "zmnc", "zmns", "numns", "numnc"):
        group.create_variable(name, (mnmax_dim,), data=getattr(surf, name))
    group.attrs["nfp"] = surf.nfp
    group.attrs["ntheta"] = surf.ntheta
    group.attrs["nzeta"] = surf.nzeta
    group.attrs["stellarator_symmetric"] = _to_bool_attr(surf.stellarator_symmetric)
    group.attrs["standard_toroidal_angle"] = _to_bool_attr(surf.standard_toroidal_angle)


def _read_fourier_surface(group):
    return dict(
        xm=np.asarray(group.variables["xm"][:]),
        xn=np.asarray(group.variables["xn"][:]),
        rmnc=np.asarray(group.variables["rmnc"][:]),
        rmns=np.asarray(group.variables["rmns"][:]),
        zmnc=np.asarray(group.variables["zmnc"][:]),
        zmns=np.asarray(group.variables["zmns"][:]),
        numns=np.asarray(group.variables["numns"][:]),
        numnc=np.asarray(group.variables["numnc"][:]),
        nfp=int(group.attrs["nfp"]),
        ntheta=int(group.attrs["ntheta"]),
        nzeta=int(group.attrs["nzeta"]),
        stellarator_symmetric=bool(group.attrs["stellarator_symmetric"]),
        standard_toroidal_angle=bool(group.attrs["standard_toroidal_angle"]),
    )


def _write_plasma(root, plasma):
    group = root.create_group("plasma")
    group.attrs["_class"] = "PlasmaSurface"
    _write_fourier_surface(group, plasma, "plasma")
    group.attrs["net_poloidal_current"] = float(plasma.net_poloidal_current)
    group.attrs["curpol"] = float(plasma.curpol)
    group.dimensions["ntheta_plasma"] = plasma.ntheta
    group.dimensions["nzeta_plasma"] = plasma.nzeta
    group.create_variable(
        "Bnormal_from_plasma_current", ("ntheta_plasma", "nzeta_plasma"),
        data=plasma.Bnormal_from_plasma_current,
    )


def _read_plasma(root):
    from .plasma_surface import PlasmaSurface

    group = root.groups["plasma"]
    data = _read_fourier_surface(group)
    plasma = PlasmaSurface(
        data["xm"], data["xn"], data["rmnc"], data["zmns"], data["rmns"], data["zmnc"],
        nfp=data["nfp"], ntheta=data["ntheta"], nzeta=data["nzeta"],
        stellarator_symmetric=data["stellarator_symmetric"],
        standard_toroidal_angle=data["standard_toroidal_angle"],
        numns=data["numns"], numnc=data["numnc"],
    )
    plasma.net_poloidal_current = float(group.attrs["net_poloidal_current"])
    plasma.curpol = float(group.attrs["curpol"])
    plasma.Bnormal_from_plasma_current = np.asarray(group.variables["Bnormal_from_plasma_current"][:])
    return plasma


def _write_coil(root, coil):
    group = root.create_group("coil")
    group.attrs["_class"] = "CoilSurface"
    _write_fourier_surface(group, coil, "coil")


def _read_coil(root):
    from .coil_surface import CoilSurface

    group = root.groups["coil"]
    data = _read_fourier_surface(group)
    return CoilSurface(
        data["xm"], data["xn"], data["rmnc"], data["zmns"], data["rmns"], data["zmnc"],
        nfp=data["nfp"], ntheta=data["ntheta"], nzeta=data["nzeta"],
        stellarator_symmetric=data["stellarator_symmetric"],
        standard_toroidal_angle=data["standard_toroidal_angle"],
        numns=data["numns"], numnc=data["numnc"],
    )


def _write_problem(root, problem):
    # Bnormal_from_net_coil_currents needs the Fortran `h` term; force the
    # (usually already-built) operators rather than saving a stale None.
    problem._ensure_operators()

    group = root.create_group("problem")
    group.attrs["_class"] = "Regcoil"
    group.attrs["mpol_potential"] = problem.mpol_potential
    group.attrs["ntor_potential"] = problem.ntor_potential
    group.attrs["net_poloidal_current"] = problem.net_poloidal_current
    group.attrs["net_toroidal_current"] = problem.net_toroidal_current
    group.attrs["stellarator_symmetric"] = _to_bool_attr(problem.stellarator_symmetric)
    group.dimensions["ntheta_plasma"] = problem.plasma.ntheta
    group.dimensions["nzeta_plasma"] = problem.plasma.nzeta
    group.create_variable(
        "Bnormal_from_net_coil_currents", ("ntheta_plasma", "nzeta_plasma"),
        data=problem.Bnormal_from_net_coil_currents,
    )


def _read_problem(root, plasma, coil):
    from .regcoil import Regcoil

    group = root.groups["problem"]
    return Regcoil._from_loaded(
        plasma, coil,
        mpol_potential=int(group.attrs["mpol_potential"]),
        ntor_potential=int(group.attrs["ntor_potential"]),
        net_poloidal_current=float(group.attrs["net_poloidal_current"]),
        net_toroidal_current=float(group.attrs["net_toroidal_current"]),
        stellarator_symmetric=bool(group.attrs["stellarator_symmetric"]),
        Bnormal_from_net_coil_currents=np.asarray(group.variables["Bnormal_from_net_coil_currents"][:]),
    )


def _write_solutions(root, solutions):
    solutions = list(solutions)
    if not solutions:
        raise ValueError("solutions is empty")
    problem = solutions[0].problem
    plasma, coil = problem.plasma, problem.coil
    nlambda = len(solutions)

    group = root.create_group("solutions")
    group.attrs["_class"] = "SolutionScan"
    group.dimensions["lambda"] = nlambda
    group.dimensions["nbf"] = problem.nbf
    group.dimensions["ntheta_plasma"] = plasma.ntheta
    group.dimensions["nzeta_plasma"] = plasma.nzeta
    group.dimensions["ntheta_coil"] = coil.ntheta
    group.dimensions["nzeta_coil"] = coil.nzeta
    group.dimensions["xyz"] = 3

    lam = np.array([sol.lam for sol in solutions])
    solution = np.stack([sol.solution for sol in solutions])
    f_B = np.array([sol.f_B for sol in solutions])
    f_K = np.array([sol.f_K for sol in solutions])
    max_K = np.array([sol.max_K for sol in solutions])
    rms_K = np.array([sol.rms_K for sol in solutions])
    max_Bnormal = np.array([sol.max_Bnormal for sol in solutions])
    Bnormal_total = np.stack([sol.Bnormal_total for sol in solutions])
    current_potential = np.stack([sol.current_potential() for sol in solutions])
    current_density = np.stack([sol.current_density() for sol in solutions])

    group.create_variable("lam", ("lambda",), data=lam)
    group.create_variable("solution", ("lambda", "nbf"), data=solution)
    group.create_variable("f_B", ("lambda",), data=f_B)
    group.create_variable("f_K", ("lambda",), data=f_K)
    group.create_variable("max_K", ("lambda",), data=max_K)
    group.create_variable("rms_K", ("lambda",), data=rms_K)
    group.create_variable("max_Bnormal", ("lambda",), data=max_Bnormal)
    group.create_variable(
        "Bnormal_total", ("lambda", "ntheta_plasma", "nzeta_plasma"), data=Bnormal_total
    )
    group.create_variable(
        "current_potential", ("lambda", "ntheta_coil", "nzeta_coil"), data=current_potential
    )
    group.create_variable(
        "current_density", ("lambda", "xyz", "ntheta_coil", "nzeta_coil"), data=current_density
    )


def _read_solutions(root, problem):
    from .regcoil import Solution, SolutionScan

    group = root.groups["solutions"]
    lam = np.asarray(group.variables["lam"][:])
    solution = np.asarray(group.variables["solution"][:])
    f_B = np.asarray(group.variables["f_B"][:])
    f_K = np.asarray(group.variables["f_K"][:])
    max_K = np.asarray(group.variables["max_K"][:])
    rms_K = np.asarray(group.variables["rms_K"][:])
    max_Bnormal = np.asarray(group.variables["max_Bnormal"][:])
    Bnormal_total = np.asarray(group.variables["Bnormal_total"][:])
    current_potential = np.asarray(group.variables["current_potential"][:])
    current_density = np.asarray(group.variables["current_density"][:])

    solutions = [
        Solution(
            problem=problem,
            lam=float(lam[i]),
            solution=solution[i],
            f_B=float(f_B[i]),
            f_K=float(f_K[i]),
            max_K=float(max_K[i]),
            rms_K=float(rms_K[i]),
            max_Bnormal=float(max_Bnormal[i]),
            Bnormal_total=Bnormal_total[i],
            _current_potential=current_potential[i],
            _current_density=current_density[i],
        )
        for i in range(lam.shape[0])
    ]
    return SolutionScan(solutions)


def save(path, *, plasma=None, coil=None, problem=None, solutions=None):
    """Write only the groups given, deduping the shared problem/surfaces:
    passing `solutions` (or `problem`) transitively carries its `problem` ->
    `coil`, `plasma`, so `save(path, solutions=scan)` writes the whole
    bundle in one call (ADR-028).
    """
    if solutions is not None:
        solutions = list(solutions)
        if not solutions:
            raise ValueError("solutions is empty")
        scan_problem = solutions[0].problem
        if any(sol.problem is not scan_problem for sol in solutions):
            raise ValueError("all solutions must share the same problem to be saved together")
        if problem is not None and problem is not scan_problem:
            raise ValueError("problem does not match solutions[0].problem")
        problem = scan_problem

    if problem is not None:
        if plasma is not None and plasma is not problem.plasma:
            raise ValueError("plasma does not match problem.plasma")
        if coil is not None and coil is not problem.coil:
            raise ValueError("coil does not match problem.coil")
        plasma = problem.plasma
        coil = problem.coil

    if plasma is None and coil is None and problem is None and solutions is None:
        raise ValueError("save() needs at least one of plasma, coil, problem, solutions")

    import h5netcdf

    with h5netcdf.File(path, "w") as f:
        f.attrs["format_version"] = FORMAT_VERSION
        from . import __version__
        f.attrs["regcoil_version"] = __version__

        if plasma is not None:
            _write_plasma(f, plasma)
        if coil is not None:
            _write_coil(f, coil)
        if problem is not None:
            _write_problem(f, problem)
        if solutions is not None:
            _write_solutions(f, solutions)


class LoadResult:
    """Container returned by `load()`. `.plasma`, `.coil`, `.problem`,
    `.solutions` are `None` if the corresponding group is absent."""

    def __init__(self, plasma=None, coil=None, problem=None, solutions=None):
        self.plasma = plasma
        self.coil = coil
        self.problem = problem
        self.solutions = solutions

    def __repr__(self):
        present = [
            name for name, value in (
                ("plasma", self.plasma), ("coil", self.coil),
                ("problem", self.problem), ("solutions", self.solutions),
            ) if value is not None
        ]
        return f"LoadResult({', '.join(present)})"


def _read_all(f):
    plasma = _read_plasma(f) if "plasma" in f.groups else None
    coil = _read_coil(f) if "coil" in f.groups else None
    problem = _read_problem(f, plasma, coil) if "problem" in f.groups else None
    solutions = _read_solutions(f, problem) if "solutions" in f.groups else None
    return plasma, coil, problem, solutions


def load(path):
    """Load whichever of `/plasma`, `/coil`, `/problem`, `/solutions` are
    present in `path`; returns a `LoadResult` with the rest `None`."""
    import h5netcdf

    with h5netcdf.File(path, "r") as f:
        plasma, coil, problem, solutions = _read_all(f)
    return LoadResult(plasma=plasma, coil=coil, problem=problem, solutions=solutions)


def load_plasma(path):
    plasma = load(path).plasma
    if plasma is None:
        raise ValueError(f"{path!r} has no /plasma group")
    return plasma


def load_coil(path):
    coil = load(path).coil
    if coil is None:
        raise ValueError(f"{path!r} has no /coil group")
    return coil


def load_problem(path):
    problem = load(path).problem
    if problem is None:
        raise ValueError(f"{path!r} has no /problem group")
    return problem


def load_solution(path):
    solutions = load(path).solutions
    if solutions is None:
        raise ValueError(f"{path!r} has no /solutions group")
    if len(solutions) != 1:
        raise ValueError(
            f"{path!r} holds {len(solutions)} solutions (a lambda scan); "
            "use regcoil.load(path).solutions to get the SolutionScan"
        )
    return solutions[0]
