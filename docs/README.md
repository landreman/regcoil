To build the docs locally, first ensure the proper packages are installed, which
can be done from the regcoil package root directory via

```bash
pip install --no-build-isolation -e .[docs]
```

Then from this directory, run

```bash
make html
```

After the compilation has completed, you can open the documentation site with

```bash
open _build/html/index.html
```
