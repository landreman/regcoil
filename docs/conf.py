# Configuration file for the Sphinx documentation builder.
# See https://www.sphinx-doc.org/en/master/usage/configuration.html

import inspect
import subprocess
import sys
from pathlib import Path

import regcoil

# -- Project information -----------------------------------------------------

project = "regcoil"
copyright = "2026, Matt Landreman"
author = "Matt Landreman"
version = regcoil.__version__
release = version

# -- General configuration ----------------------------------------------------

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.autosummary",
    "sphinx.ext.napoleon",
    "sphinx.ext.doctest",
    "sphinx.ext.linkcode",
    "sphinx.ext.mathjax",
    "sphinx_copybutton",
    "myst_nb",
]

templates_path = ["_templates"]
exclude_patterns = [
    "_build",
    "Thumbs.db",
    ".DS_Store",
    ".jupyter_cache",
    "jupyter_execute",
    # Internal migration/engineering planning docs, not part of the public
    # user-facing manual (see docs/migration/DECISIONS.md ADR-030).
    "migration",
    # Developer-facing build instructions (this README).
    "README.md",
]

# MyST / notebook sources: both plain prose pages and executed tutorial pages
# are authored as MyST Markdown (jupytext-compatible; see ADR-030). Any page
# containing a ``{code-cell}`` directive is treated as a notebook and executed.
# (Do not set `source_suffix` explicitly: myst_nb registers the `.md` ->
# MyST-NB source parser itself; overriding it here maps `.md` back to plain
# reStructuredText and silently breaks every `{directive}` fence.)
myst_enable_extensions = ["dollarmath", "colon_fence"]
myst_heading_anchors = 3

# -- MyST-NB (executed docs, ADR-030) -----------------------------------------

nb_execution_mode = "cache"
nb_execution_raise_on_error = True
nb_execution_timeout = 300
nb_execution_show_tb = True
# Plain prose pages (no code-cell) are not notebooks; nothing to execute there.

# -- Autosummary / autodoc / napoleon ------------------------------------------

autosummary_generate = True
autosummary_imported_members = False
autodoc_typehints = "description"
napoleon_numpy_docstring = True
napoleon_google_docstring = False
napoleon_use_rtype = False

# -- doctest (README quickstart snippet, ADR-030) ------------------------------

doctest_global_setup = "import regcoil"

# -- linkcode: per-object [source] links to GitHub -----------------------------

_GITHUB_URL = "https://github.com/landreman/pyREGCOIL"
_REPO_ROOT = Path(__file__).resolve().parent.parent


def _git_head() -> str:
    try:
        return (
            subprocess.check_output(["git", "log", "-n1", "--pretty=%H"], cwd=_REPO_ROOT)
            .strip()
            .decode("utf-8")
        )
    except subprocess.CalledProcessError:
        return "master"


_BLOB = _git_head()


def linkcode_resolve(domain, info):
    """Return a GitHub URL for a documented Python object, or None."""
    if domain != "py" or not info["module"]:
        return None

    modname = info["module"]
    fullname = info["fullname"]

    submod = sys.modules.get(modname)
    if submod is None:
        return None

    obj = submod
    for part in fullname.split("."):
        try:
            obj = getattr(obj, part)
        except AttributeError:
            return None

    obj = inspect.unwrap(obj)

    try:
        filename = inspect.getsourcefile(obj)
        source, lineno = inspect.getsourcelines(obj)
    except (TypeError, OSError):
        # e.g. objects implemented in the compiled `regcoil._core` extension.
        return None

    if filename is None:
        return None

    try:
        rel_path = Path(filename).resolve().relative_to(_REPO_ROOT)
    except ValueError:
        return None

    start = lineno
    end = lineno + len(source) - 1
    return f"{_GITHUB_URL}/blob/{_BLOB}/{rel_path.as_posix()}#L{start}-L{end}"


# -- HTML output ----------------------------------------------------------------

html_theme = "furo"
html_title = "regcoil"
