"""Unit tests for `Surface.scale` / `PlasmaSurface.scale`."""

from pathlib import Path

import numpy as np
import pytest

from regcoil import CoilSurface, FourierSurface, PlasmaSurface, Regcoil

REPO_ROOT = Path(__file__).resolve().parents[2]
EQUILIBRIA = REPO_ROOT / "equilibria"
WOUT = EQUILIBRIA / "wout_li383_1.4m.nc"


def _torus(cls=FourierSurface, R0=3.0, a=0.7, nfp=3, ntheta=32, nzeta=30):
    return cls.circular_torus(R0=R0, a=a, nfp=nfp, ntheta=ntheta, nzeta=nzeta)


def test_length_factor_scales_geometry():
    surf = _torus()
    # Touch the caches first, so a stale one would be caught below.
    area, volume, aspect = surf.area, surf.volume, surf.aspect_ratio
    r_before = surf.r.copy()

    scaled = surf.scale(length_factor=2.0)

    np.testing.assert_allclose(scaled.r, 2.0 * r_before, rtol=1e-13)
    np.testing.assert_allclose(scaled.minor_radius, 2 * surf.minor_radius, rtol=1e-13)
    np.testing.assert_allclose(scaled.major_radius, 2 * surf.major_radius, rtol=1e-13)
    np.testing.assert_allclose(scaled.area, 4 * area, rtol=1e-13)
    np.testing.assert_allclose(scaled.volume, 8 * volume, rtol=1e-13)
    np.testing.assert_allclose(scaled.aspect_ratio, aspect, rtol=1e-13)

    # The original is untouched.
    np.testing.assert_allclose(surf.r, r_before, rtol=0, atol=0)
    np.testing.assert_allclose(surf.area, area, rtol=0, atol=0)


@pytest.mark.parametrize("target", ["minor_radius", "major_radius", "volume", "area"])
def test_length_targets_are_hit(target):
    surf = _torus()
    wanted = 1.7 * getattr(surf, target)
    scaled = surf.scale(**{target: wanted})
    np.testing.assert_allclose(getattr(scaled, target), wanted, rtol=1e-12)
    # Only the size changed, not the shape.
    np.testing.assert_allclose(scaled.aspect_ratio, surf.aspect_ratio, rtol=1e-12)


def test_no_target_is_identity_but_a_copy():
    surf = _torus()
    same = surf.scale()
    assert same is not surf
    np.testing.assert_allclose(same.r, surf.r, rtol=0, atol=0)


def test_subclass_and_attributes_are_preserved():
    coil = _torus(CoilSurface)
    scaled = coil.scale(length_factor=1.5)
    assert type(scaled) is CoilSurface
    assert scaled.nfp == coil.nfp
    assert scaled.ntheta == coil.ntheta
    assert scaled.nzeta == coil.nzeta
    assert scaled.stellarator_symmetric == coil.stellarator_symmetric


def test_nu_modes_are_not_scaled():
    """`nu` is an angle, so a surface whose `zeta` is not the standard toroidal
    angle must scale to exactly `factor * r` -- which it only does if `numns`
    is left alone."""
    plasma = PlasmaSurface.from_wout(str(WOUT), ntheta=16, nzeta=15)
    coil = CoilSurface.from_uniform_offset(
        plasma, separation=0.3, ntheta=16, nzeta=15, mpol=6, ntor=5,
        standard_toroidal_angle=False,
    )
    assert not coil.standard_toroidal_angle
    assert np.any(coil.numns != 0)

    scaled = coil.scale(length_factor=3.0)
    np.testing.assert_allclose(scaled.numns, coil.numns, rtol=0, atol=0)
    np.testing.assert_allclose(scaled.numnc, coil.numnc, rtol=0, atol=0)
    np.testing.assert_allclose(scaled.r, 3.0 * coil.r, rtol=1e-13)


def test_plasma_length_scaling():
    plasma = PlasmaSurface.from_wout(str(WOUT), ntheta=16, nzeta=15)
    c = 2.5
    scaled = plasma.scale(length_factor=c)

    # ~ B*L, so one power of the length factor:
    np.testing.assert_allclose(
        scaled.net_poloidal_current, c * plasma.net_poloidal_current, rtol=1e-13
    )
    np.testing.assert_allclose(scaled.curpol, c * plasma.curpol, rtol=1e-13)
    # The fields themselves are length-invariant:
    np.testing.assert_allclose(scaled.modB, plasma.modB, rtol=0, atol=0)
    np.testing.assert_allclose(
        scaled.volume_averaged_B, plasma.volume_averaged_B, rtol=0, atol=0
    )
    np.testing.assert_allclose(
        scaled.Bnormal_from_plasma_current,
        plasma.Bnormal_from_plasma_current,
        rtol=0, atol=0,
    )


def test_plasma_field_scaling():
    plasma = PlasmaSurface.from_wout(str(WOUT), ntheta=16, nzeta=15)
    plasma.set_bnormal_from_bnorm_file(str(EQUILIBRIA / "bnorm.li383_1.4m"))
    assert np.any(plasma.Bnormal_from_plasma_current != 0)
    r_before = plasma.r.copy()
    c = 1.75

    scaled = plasma.scale(B_factor=c)

    np.testing.assert_allclose(
        scaled.net_poloidal_current, c * plasma.net_poloidal_current, rtol=1e-13
    )
    np.testing.assert_allclose(scaled.curpol, c * plasma.curpol, rtol=1e-13)
    np.testing.assert_allclose(scaled.modB, c * plasma.modB, rtol=1e-13)
    np.testing.assert_allclose(
        scaled.volume_averaged_B, c * plasma.volume_averaged_B, rtol=1e-13
    )
    np.testing.assert_allclose(
        scaled.Bnormal_from_plasma_current,
        c * plasma.Bnormal_from_plasma_current,
        rtol=1e-13,
    )
    # Geometry untouched.
    np.testing.assert_allclose(scaled.r, r_before, rtol=0, atol=0)
    np.testing.assert_allclose(plasma.r, r_before, rtol=0, atol=0)


def test_plasma_combined_targets():
    plasma = PlasmaSurface.from_wout(str(WOUT), ntheta=16, nzeta=15)
    a_new, B_new = 1.7, 5.7
    cL = a_new / plasma.minor_radius
    cB = B_new / plasma.volume_averaged_B

    scaled = plasma.scale(minor_radius=a_new, volume_averaged_B=B_new)

    np.testing.assert_allclose(scaled.minor_radius, a_new, rtol=1e-12)
    np.testing.assert_allclose(scaled.volume_averaged_B, B_new, rtol=1e-13)
    np.testing.assert_allclose(
        scaled.net_poloidal_current, cL * cB * plasma.net_poloidal_current, rtol=1e-12
    )
    np.testing.assert_allclose(scaled.curpol, cL * cB * plasma.curpol, rtol=1e-12)
    np.testing.assert_allclose(scaled.modB, cB * plasma.modB, rtol=1e-13)


def test_max_modB_target():
    plasma = PlasmaSurface.from_wout(str(WOUT), ntheta=16, nzeta=15)
    wanted = 3.3
    scaled = plasma.scale(max_modB=wanted)
    np.testing.assert_allclose(scaled.modB.max(), wanted, rtol=1e-13)


def test_scale_length_and_scale_B_shorthands():
    plasma = PlasmaSurface.from_wout(str(WOUT), ntheta=16, nzeta=15)
    np.testing.assert_allclose(
        plasma.scale_length(2.0).r, plasma.scale(length_factor=2.0).r, rtol=0, atol=0
    )
    np.testing.assert_allclose(
        plasma.scale_B(2.0).modB, plasma.scale(B_factor=2.0).modB, rtol=0, atol=0
    )


def test_bnormal_from_coil_currents_scales_as_the_field_factor():
    """The physics check: scale the plasma by (cL, cB) and the coil surface by
    cL, and the field the net coil currents produce on the plasma boundary
    must come out cB times larger -- the whole point of the dimensional
    bookkeeping in `net_poloidal_current`. (A shaped, non-axisymmetric pair is
    needed here: for circular tori this field is zero to roundoff.)
    """
    cL, cB = 2.0, 1.5
    plasma = PlasmaSurface.from_wout(str(WOUT), ntheta=12, nzeta=11)
    coil = CoilSurface.from_nescin(
        str(EQUILIBRIA / "nescin.li383_realWindingSurface"), nfp=plasma.nfp,
        ntheta=12, nzeta=11,
    )

    def build(p, c):
        return Regcoil(
            p, c, mpol_potential=3, ntor_potential=2,
            net_toroidal_current=0.0, stellarator_symmetric=True,
        )

    base = build(plasma, coil)
    scaled = build(
        plasma.scale(length_factor=cL, B_factor=cB), coil.scale(length_factor=cL)
    )

    scale_of_B = np.abs(base.Bnormal_from_net_coil_currents).max()
    assert scale_of_B > 1e-3
    np.testing.assert_allclose(
        scaled.Bnormal_from_net_coil_currents,
        cB * base.Bnormal_from_net_coil_currents,
        rtol=1e-11, atol=1e-13 * scale_of_B,
    )


def test_scaled_plasma_survives_a_save_load_roundtrip(tmp_path):
    plasma = PlasmaSurface.from_wout(str(WOUT), ntheta=16, nzeta=15).scale(
        minor_radius=1.7, volume_averaged_B=5.7
    )
    path = tmp_path / "scaled.nc"
    plasma.save(str(path))
    loaded = PlasmaSurface.load(str(path))
    np.testing.assert_allclose(loaded.r, plasma.r, rtol=1e-13)
    np.testing.assert_allclose(
        loaded.net_poloidal_current, plasma.net_poloidal_current, rtol=1e-13
    )
    np.testing.assert_allclose(loaded.curpol, plasma.curpol, rtol=1e-13)
    np.testing.assert_allclose(
        loaded.volume_averaged_B, plasma.volume_averaged_B, rtol=1e-13
    )


def test_errors():
    surf = _torus()
    plasma = PlasmaSurface.circular_torus(R0=3.0, a=0.7, nfp=3, ntheta=16, nzeta=15)

    with pytest.raises(ValueError, match="At most one length target"):
        surf.scale(minor_radius=1.0, volume=1.0)
    with pytest.raises(ValueError, match="At most one field target"):
        plasma.scale(B_factor=2.0, max_modB=1.0)
    with pytest.raises(ValueError, match="must be positive"):
        surf.scale(length_factor=-1.0)
    with pytest.raises(ValueError, match="must be positive"):
        surf.scale(volume=0.0)
    with pytest.raises(ValueError, match="no volume_averaged_B"):
        plasma.scale(volume_averaged_B=5.0)
    with pytest.raises(ValueError, match="modB is identically zero"):
        plasma.scale(max_modB=5.0)
    # Field targets are not available on a surface with no field.
    with pytest.raises(TypeError):
        surf.scale(B_factor=2.0)
