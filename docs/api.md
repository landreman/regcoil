# API reference

## Surfaces

```{eval-rst}
.. autosummary::
   :toctree: _api/
   :recursive:

   regcoil.Surface
   regcoil.FourierSurface
   regcoil.PlasmaSurface
   regcoil.CoilSurface
```

## Theta reparameterization

Reparameterize a surface's poloidal angle without changing its shape --
`Surface.reparameterize_theta`, or
`CoilSurface.from_uniform_offset(theta_reparameterization=...)`.

```{eval-rst}
.. autosummary::
   :toctree: _api/
   :recursive:

   regcoil.UniformArclength
   regcoil.CurvatureWeighted
   regcoil.ThetaMap
```

## Problem and solution

```{eval-rst}
.. autosummary::
   :toctree: _api/
   :recursive:

   regcoil.Regcoil
   regcoil.Solution
   regcoil.SolutionScan
```

## Plotting

```{eval-rst}
.. autosummary::
   :toctree: _api/
   :recursive:

   regcoil.plot.cross_sections
   regcoil.plot.cross_sections_overlay
   regcoil.plot.pareto
   regcoil.plot.lambda_scan
   regcoil.plot.current_potential
   regcoil.plot.current_potential_scan
   regcoil.plot.current_density
   regcoil.plot.bnormal
   regcoil.plot.bnormal_scan
   regcoil.plot.plot_3d
   regcoil.plot.coil_3d
   regcoil.plot.all
```

## Magnetic fields and field lines

```{eval-rst}
.. autosummary::
   :toctree: _api/
   :recursive:

   regcoil.MagneticFieldEvaluator
   regcoil.fieldlines.ToroidalTracer
   regcoil.fieldlines.ToroidalTracingOptions
   regcoil.fieldlines.ToroidalDiagnostic
   regcoil.fieldlines.PoincareTrace
   regcoil.fieldlines.PoincareBatch
   regcoil.fieldlines.MagneticAxisResult
   regcoil.fieldlines.IotaProfile
   regcoil.fieldlines.cylindrical_field
   regcoil.fieldlines.estimate_iota
   regcoil.plot.poincare
   regcoil.plot.iota_profile
```

## Cutting into discrete coils

```{eval-rst}
.. autosummary::
   :toctree: _api/
   :recursive:

   regcoil.cut.cut
   regcoil.CutCoils
```

## Saving and loading

```{eval-rst}
.. autosummary::
   :toctree: _api/
   :recursive:

   regcoil.save
   regcoil.load
```

## Bundled example datasets

```{eval-rst}
.. autosummary::
   :toctree: _api/
   :recursive:

   regcoil.examples
   regcoil.ExampleDataset
```

## Utilities

```{eval-rst}
.. autosummary::
   :toctree: _api/
   :recursive:

   regcoil.log
   regcoil.omp_max_threads
```
