"""Fast synthetic tests for the public toroidal field-line API."""

import numpy as np
import pytest

from regcoil.fieldlines import (
    ToroidalTracer,
    ToroidalTracingOptions,
    cylindrical_field,
    estimate_iota,
)


AXIS = np.array([6.0, 0.15])
IOTA = 0.37


def rotating_field(*, sign=1.0, iota=IOTA):
    """Axisymmetric field whose cylindrical field-line RHS is a rotation."""

    def field(point):
        x, y, z = np.asarray(point)
        radius = np.hypot(x, y)
        cosine = x / radius
        sine = y / radius
        B_phi = float(sign)
        B_R = B_phi * (-iota * (z - AXIS[1])) / radius
        B_Z = B_phi * (iota * (radius - AXIS[0])) / radius
        return np.array([
            cosine * B_R - sine * B_phi,
            sine * B_R + cosine * B_phi,
            B_Z,
        ])

    return field


class CircularPlasma:
    def cross_section(self, phi=0.0):
        theta = np.linspace(0.0, 2 * np.pi, 129)
        return AXIS[0] + np.cos(theta), AXIS[1] + np.sin(theta)


def test_cylindrical_field_round_trip():
    phi = 0.41
    rz = np.array([6.2, -0.1])
    expected = rotating_field()(np.array([
        rz[0] * np.cos(phi), rz[0] * np.sin(phi), rz[1]
    ]))
    B_R, B_phi, B_Z = cylindrical_field(phi, rz, rotating_field())
    np.testing.assert_allclose(
        [B_R, B_phi, B_Z],
        [
            np.cos(phi) * expected[0] + np.sin(phi) * expected[1],
            -np.sin(phi) * expected[0] + np.cos(phi) * expected[1],
            expected[2],
        ],
        atol=2e-15,
    )


def test_successful_diagnostic_evaluates_magnetic_field_once():
    evaluations = 0
    field = rotating_field()

    def counting_field(point):
        nonlocal evaluations
        evaluations += 1
        return field(point)

    diagnostic = ToroidalTracer(counting_field, nfp=1).diagnose(
        0.23, [6.2, -0.1]
    )
    assert diagnostic.valid
    assert evaluations == 1


@pytest.mark.parametrize("sign", [1.0, -1.0])
def test_sections_and_iota_work_in_both_toroidal_directions(sign):
    tracer = ToroidalTracer(
        rotating_field(sign=sign),
        nfp=3,
        options=ToroidalTracingOptions(rtol=2e-11, atol=2e-13),
    )
    initial = np.array([[6.15, AXIS[1]], [6.3, AXIS[1]]])
    phases = np.array([0.0, 0.13, 0.5, 0.82])
    batch = tracer.trace_poincare(
        initial,
        periods=40,
        phase_fractions=phases,
        direction="auto",
        max_workers=2,
    )

    assert batch.direction == int(sign)
    assert len(batch.successful) == 2
    assert [trace.initial_rz[0] for trace in batch.traces] == [6.15, 6.3]
    for trace in batch.successful:
        for phase_index, expected_phase in enumerate(phases):
            sample_phi, _ = trace.section(phase_index)
            np.testing.assert_allclose(
                np.mod(sample_phi / tracer.field_period, 1.0),
                expected_phase,
                atol=2e-14,
            )
        return_phi, return_rz = trace.same_section_returns()
        assert len(return_phi) == 41
        assert return_rz.shape == (41, 2)

    profile = estimate_iota(batch, AXIS, period_counts=(10, 20, 40))
    np.testing.assert_allclose(profile.iota, IOTA, rtol=2e-9, atol=2e-9)
    np.testing.assert_allclose(profile.radius, [0.15, 0.3])


def test_expected_guard_failure_is_a_record_not_a_batch_exception():
    def field_without_toroidal_component(point):
        return np.array([1.0, 0.0, 0.0])

    tracer = ToroidalTracer(field_without_toroidal_component, nfp=1)
    batch = tracer.trace_poincare([[6.0, 0.0]], periods=2)
    assert not batch.traces[0].success
    assert "absolute floor" in batch.traces[0].message
    assert batch.traces[0].failure_diagnostic is not None


def test_auto_direction_rejects_mixed_initial_bphi_signs():
    def mixed_field(point):
        sign = 1.0 if point[0] < 6.0 else -1.0
        return np.array([0.0, sign, 0.0])

    tracer = ToroidalTracer(mixed_field, nfp=1)
    with pytest.raises(ValueError, match="mixed signs"):
        tracer.trace_poincare(
            [[5.5, 0.0], [6.5, 0.0]], periods=1, direction="auto"
        )


def test_finds_elliptic_axis_of_synthetic_return_map():
    tracer = ToroidalTracer(rotating_field(), nfp=1)
    axis = tracer.find_magnetic_axis(CircularPlasma())
    np.testing.assert_allclose(axis.rz, AXIS, atol=2e-9)
    assert axis.return_error < 1e-8
    assert axis.determinant == pytest.approx(1.0, abs=2e-6)
    assert np.allclose(np.abs(axis.eigenvalues), 1.0, atol=2e-6)


def test_axis_options_inherit_tracer_guards_and_tighten_tolerances(monkeypatch):
    tracer_options = ToroidalTracingOptions(
        bphi_absolute_floor=2.0e-8,
        bphi_relative_floor=2.0e-2,
        maximum_rhs_norm=23.0,
        minimum_field_magnitude=4.0e-12,
        method="DOP853",
        rtol=3.0e-6,
        atol=4.0e-8,
    )
    tracer = ToroidalTracer(
        rotating_field(), nfp=1, options=tracer_options
    )
    observed_options = []
    angle = 0.4
    return_map_jacobian = np.array([
        [np.cos(angle), -np.sin(angle)],
        [np.sin(angle), np.cos(angle)],
    ])

    def fake_first_return_map(rz, *, direction=1, options=None):
        observed_options.append(options)
        return AXIS + return_map_jacobian @ (np.asarray(rz) - AXIS)

    monkeypatch.setattr(tracer, "first_return_map", fake_first_return_map)
    tracer.find_magnetic_axis(CircularPlasma())

    assert observed_options
    axis_options = observed_options[0]
    assert all(options is axis_options for options in observed_options)
    assert axis_options.bphi_absolute_floor == 2.0e-8
    assert axis_options.bphi_relative_floor == 2.0e-2
    assert axis_options.maximum_rhs_norm == 23.0
    assert axis_options.minimum_field_magnitude == 4.0e-12
    assert axis_options.method == "DOP853"
    assert axis_options.rtol == 1.0e-10
    assert axis_options.atol == 1.0e-11


@pytest.mark.parametrize(
    "kwargs, message",
    [
        ({"periods": 0}, "periods"),
        ({"periods": 1, "phase_fractions": (0.25,)}, "contain 0"),
        ({"periods": 1, "phase_fractions": (0.0, 1.0)}, "distinct"),
        ({"periods": 1, "max_workers": 0}, "max_workers"),
    ],
)
def test_rejects_invalid_batch_options(kwargs, message):
    tracer = ToroidalTracer(rotating_field(), nfp=1)
    with pytest.raises(ValueError, match=message):
        tracer.trace_poincare([[6.2, AXIS[1]]], **kwargs)
