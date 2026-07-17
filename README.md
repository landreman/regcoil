### REGCOIL

A regularized current potential method for fast computation of the shapes of stellarator coils

![](https://github.com/landreman/regcoil/blob/master/manual/m20170111_01_compareNescoilToRegcoilCoils.png)

This program is described in the paper 
`M Landreman, "An improved current potential method for fast computation of stellarator coil shapes," Nuclear Fusion 57, 046003 (2017)`,
which is available in this repository and also [at arXiv:1609.04378](https://arxiv.org/pdf/1609.04378.pdf).
For further documentation, see the user manual [here](http://landreman.github.io/regcoil/regcoilManual.pdf).

Python packaging / migration plans: [`docs/migration/`](docs/migration/OVERVIEW.md).

### Build and test (legacy executable)

```bash
make
make test
```

On macOS with Homebrew NetCDF (Apple Silicon), the default makefile paths under `/usr/local` are often wrong—see [`docs/migration/LOCAL.md`](docs/migration/LOCAL.md) for the `EXTRA_*` flag overrides that work locally.
