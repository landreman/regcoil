"""Smoke test: package is importable after install / editable layout."""

from io import StringIO
import re

# packaging arrives transitively with matplotlib, a hard runtime dependency.
from packaging.version import Version

import regcoil


def test_log_enables_timestamped_package_logging():
    stream = StringIO()
    regcoil.log(stream=stream)

    try:
        import logging

        logging.getLogger("regcoil.regcoil").info("hello")
        output = stream.getvalue().strip()
        lines = output.splitlines()
        assert re.fullmatch(r"\d{2}:\d{2}:\d{2}\.\d{2} REGCOIL: starting package logging", lines[0])
        assert re.fullmatch(r"\d{2}:\d{2}:\d{2}\.\d{2} OpenMP max threads: \d+", lines[1])
        assert re.fullmatch(r"\d{2}:\d{2}:\d{2}\.\d{2} hello", lines[2])
    finally:
        package_logger = logging.getLogger("regcoil")
        package_logger.handlers.clear()
        package_logger.propagate = True


def test_omp_max_threads():
    assert regcoil.omp_max_threads() >= 1


def test_version_matches_installed_distribution_metadata():
    """regcoil.__version__ must agree with what pip/PyPI see.

    The version literal lives only in pyproject.toml; __init__ reads it back
    from the installed distribution.  This guards against a build that somehow
    ships mismatched metadata, and documents the invariant the release workflow
    also checks against the git tag.
    """
    from importlib.metadata import version

    assert regcoil.__version__ == version("regcoil")
    # A PEP 440 parse is the same check PyPI will apply at upload time.
    Version(regcoil.__version__)
