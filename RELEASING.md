# Releasing regcoil

Maintainer notes. Users do not need this file.

## How the version is wired

There is exactly **one** version literal in the repository:

```toml
# pyproject.toml
[project]
version = "0.1.0"
```

Everything else derives from it:

| Consumer | How it gets the version |
|---|---|
| wheel / sdist metadata, PyPI | meson-python reads `[project] version` directly |
| `regcoil.__version__` | `importlib.metadata.version("regcoil")`, i.e. read back out of the installed distribution |
| `docs/conf.py` | `regcoil.__version__` |
| HDF5 output `regcoil_version` attribute | `regcoil.__version__` |
| the git tag | written by `tbump`, then re-checked by CI at release time |

`meson.build`'s `project(version: '0.0.0')` is a **placeholder that is never used** —
meson-python only consults it when `pyproject.toml` declares the version dynamic,
which it does not. Leave it at `0.0.0`.

Two independent guards keep this honest:

- `tests/unit/test_import.py::test_version_matches_installed_distribution_metadata`
  fails if the import-time version ever disagrees with the installed metadata.
- The `check-version` job in `.github/workflows/release.yml` refuses to build or
  publish anything if the release tag does not match `pyproject.toml`.

## One-time setup (manual — see the checklist at the bottom)

Publishing uses **PyPI Trusted Publishing (OIDC)**. There is no API token stored
in GitHub secrets and nothing to rotate.

## Cutting a release

### 1. Rehearse against TestPyPI

From the GitHub Actions tab, run the **Release** workflow via *Run workflow*
(`workflow_dispatch`) with *Upload to TestPyPI* checked. This builds the full
wheel matrix plus the sdist, smoke-tests every wheel, and uploads to TestPyPI —
without creating a tag or touching real PyPI.

Then install the result on a clean machine and confirm it works:

```bash
python -m venv /tmp/rc && /tmp/rc/bin/pip install -U pip
/tmp/rc/bin/pip install --index-url https://test.pypi.org/simple/ \
    --extra-index-url https://pypi.org/simple/ regcoil
/tmp/rc/bin/python -c "import regcoil; print(regcoil.__version__); print(regcoil.omp_max_threads())"
```

The `--extra-index-url` is required: TestPyPI does not mirror numpy/scipy/etc.

For the very first release, do this with a pre-release version (`tbump 0.1.0rc1`)
so a botched rehearsal does not burn the `0.1.0` version number — **a version can
never be re-uploaded to PyPI, even after deletion.**

### 2. Bump and tag

```bash
git checkout master && git pull      # tbump pushes to the current branch's upstream
tbump 0.2.0
```

`tbump` will show you exactly what it plans to do and ask for confirmation. It
rewrites the version in `pyproject.toml`, runs the fast test suite as a
pre-commit hook, commits as `Release 0.2.0`, creates the annotated tag `v0.2.0`,
and pushes the branch and tag together.

Run it from a branch that has an upstream (normally `master`). On a local-only
branch, tbump renders an incomplete `git push` command and the push fails.

### 3. Publish the GitHub Release

Go to *Releases → Draft a new release*, choose the existing `v0.2.0` tag, use
*Generate release notes*, edit as needed, and **Publish**.

Publishing is what triggers the build-and-upload. Within ~20 minutes the wheels
and sdist appear on PyPI.

### 4. Verify

```bash
pip download --no-deps --only-binary :all: regcoil==0.2.0
```

## Failure modes and what to do

**The `check-version` job fails.** The tag and `pyproject.toml` disagree. Delete
the GitHub release *and* the tag (`git push --delete origin v0.2.0`), fix the
version with `tbump`, and re-release. Nothing was uploaded, so the version number
is still usable.

**A wheel build fails after the release is published.** Nothing has been uploaded
(the publish job runs only after the whole matrix succeeds), but the GitHub
release exists. Fix the problem, delete the release and tag, and re-cut. Since
nothing reached PyPI, you may reuse the version number.

**The upload partially succeeded.** PyPI uploads are effectively atomic per
file, but if you end up with some files uploaded and some not, you cannot
overwrite them — bump to the next patch version and release again.

## What is deliberately *not* built

- **Windows wheels.** Mixing MinGW `gfortran` with MSVC-built CPython is a
  substantial infrastructure project. Windows users fall back to the sdist, or
  better, to conda. Revisit if there is real demand.
- **musllinux (Alpine) wheels.** Would need a separate apk toolchain for
  gfortran and BLAS.
- **PyPy wheels.** The extension is written against the CPython/NumPy C API.
- **Free-threaded (`t`) builds.** Not exercised; the extension uses OpenMP.

The sdist covers all of these for anyone with a Fortran compiler and BLAS.

## Platform notes for the wheel builds

- **Linux** uses the `manylinux_2_28` (AlmaLinux 8) images, with
  `gcc-gfortran` and `openblas-devel` installed via `dnf` in `before-all`.
  `auditwheel` vendors `libgfortran`, `libgomp` and `libopenblas` into the wheel.
  arm64 builds run on native `ubuntu-24.04-arm` runners — no QEMU emulation.
- **macOS** uses Homebrew's gcc, matching the CI workflow, and `delocate` vendors
  `libgfortran` and `libgomp`. BLAS comes from the system Accelerate framework
  and needs no vendoring.

  `MACOSX_DEPLOYMENT_TARGET` is set per-runner in `release.yml` (14.0 for arm64
  on `macos-14`, 15.0 for x86_64 on `macos-15-intel`). **It must be at least the
  SDK version Homebrew's gfortran was
  built against**, or delocate will refuse the wheel. This is the reason the
  wheels do not support older macOS; lowering the floor would mean switching to
  a gfortran built against an older SDK (as SciPy does).
- **OpenMP caveat:** the wheels vendor `libgomp`. If a user's NumPy/SciPy brings
  a different OpenMP runtime into the same process, oversubscription or (rarely)
  a crash is possible. This is the standard situation for OpenMP-using wheels and
  has not been observed in practice here.

---

## Manual setup checklist

These require account access and cannot be scripted from the repository.

- [ ] **Reserve the name on PyPI.** `regcoil` was unclaimed as of this writing.
      Trusted Publishing can create the project on first upload (see below), so
      no placeholder upload is needed — but the name is first-come, first-served.

- [ ] **Create the PyPI trusted publisher.** On PyPI → *Your projects* → (or
      *Publishing* → *Add a pending publisher* if the project does not exist yet):
      - PyPI Project Name: `regcoil`
      - Owner: `landreman`
      - Repository name: `regcoil`
      - Workflow name: `release.yml`
      - Environment name: `pypi`

- [ ] **Create the TestPyPI trusted publisher.** Same values, on
      test.pypi.org, but environment name `testpypi`.

- [ ] **Create the two GitHub environments.** Repo *Settings → Environments*:
      `pypi` and `testpypi`. The names must match the workflow exactly.
      Recommended: on the `pypi` environment, add yourself as a *required
      reviewer*. That gives you a final human confirmation between "wheels built
      successfully" and "uploaded to PyPI forever."

- [ ] **First release only:** rehearse with `0.1.0rc1` end-to-end before cutting
      `0.1.0`.
