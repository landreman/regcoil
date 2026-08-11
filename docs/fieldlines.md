# Field-line tracing

REGCOIL can prepare the continuous winding-surface-current magnetic field and
trace it using physical toroidal angle. This field does **not** include a
plasma-current contribution or the field of discrete cut coils.

```python
import numpy as np
import regcoil

data = regcoil.load("regcoil_out.example.nc")
solution = data.solutions[0]
tracer = regcoil.fieldlines.ToroidalTracer.from_solution(solution)

initial_R = np.linspace(5.9, 6.3, 24)
batch = tracer.trace_poincare(
    np.column_stack([initial_R, np.zeros_like(initial_R)]),
    periods=200,
    phase_fractions=(0, 0.25, 0.5, 0.75),
    direction="auto",
    max_workers=3,
)
fig = regcoil.plot.poincare(
    batch, plasma=data.plasma, winding_surface=data.coil
)
```

`ToroidalTracer.from_solution` prepares and caches the full-torus source data
once. Every adaptive ODE stage then uses the fused singleton field kernel.
Uniform `theta_stride` and `zeta_stride` values greater than one can accelerate
the field evaluation, but they require a convergence study for the intended
field-line diagnostic.

Toroidal-angle tracing is rejected when $B_\phi$ is too small, when
$|B_\phi|/|B|$ indicates poor conditioning, or when the cylindrical RHS is too
large. Expected failures are retained in `batch.failed`; unrelated programming
errors propagate normally.

An elliptic fixed point and a geometric rotational-transform profile can be
computed from the same tracer:

```python
axis = tracer.find_magnetic_axis(data.plasma, direction=batch.direction)
profile = regcoil.fieldlines.estimate_iota(
    batch, axis.rz, period_counts=(24, 50, 100, 200)
)
ax = regcoil.plot.iota_profile(profile)
```

The iota estimator deliberately uses repeated returns to the physical
$\phi=0$ section. Samples from other phases have a different magnetic-axis
position and cannot share the same geometric poloidal-angle reference. The
reported radius is geometric distance from the fitted axis, not normalized
toroidal flux.
