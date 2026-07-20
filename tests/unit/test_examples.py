"""Unit tests for `regcoil.examples`: the registry of bundled example datasets
(VMEC `wout`, BNORM, and nescin winding-surface files installed into the package
as `regcoil/equilibria/`). See `src/regcoil/examples.py`.
"""

import dataclasses
from pathlib import Path

import pytest

import regcoil
from regcoil import ExampleDataset, examples


# canonical id -> the exact filenames the dataset must resolve to
_EXPECTED_FILES = {
    "li383_1.4m": (
        "wout_li383_1.4m.nc",
        "bnorm.li383_1.4m",
        "nescin.li383_realWindingSurface",
    ),
    "d23p4_tm": (
        "wout_d23p4_tm.nc",
        "bnorm.d23p4_tm",
        "nescin.w7x_winding_surface_from_Drevlak",
    ),
}


def test_exported_from_top_level_package():
    assert regcoil.examples is examples
    assert regcoil.ExampleDataset is ExampleDataset
    assert "examples" in regcoil.__all__
    assert "ExampleDataset" in regcoil.__all__


def test_available_lists_names_and_aliases():
    names = examples.available()
    assert names == ["NCSX", "li383", "li383_1.4m", "W7-X", "d23p4_tm"]
    # available() returns a fresh copy each call; mutating it must not stick.
    names.append("bogus")
    assert "bogus" not in examples.available()


@pytest.mark.parametrize(
    "alias, canonical",
    [
        ("NCSX", "li383_1.4m"),
        ("li383", "li383_1.4m"),
        ("li383_1.4m", "li383_1.4m"),
        ("W7-X", "d23p4_tm"),
        ("W7X", "d23p4_tm"),
        ("d23p4_tm", "d23p4_tm"),
    ],
)
def test_aliases_resolve_to_canonical(alias, canonical):
    assert examples(alias).name == canonical


@pytest.mark.parametrize("raw", ["ncsx", "  NcSx  ", "\tw7-x\n"])
def test_lookup_is_case_insensitive_and_stripped(raw):
    # These normalize to a known dataset rather than raising.
    assert examples(raw).name in _EXPECTED_FILES


@pytest.mark.parametrize("canonical", list(_EXPECTED_FILES))
def test_resolved_files_exist_and_have_expected_names(canonical):
    ds = examples(canonical)
    wout_name, bnorm_name, nescin_name = _EXPECTED_FILES[canonical]

    assert ds.wout.name == wout_name
    assert ds.bnorm.name == bnorm_name
    assert ds.nescin.name == nescin_name

    for path in (ds.wout, ds.bnorm, ds.nescin):
        assert isinstance(path, Path)
        assert path.is_absolute()
        assert path.is_file()
        # Resolved via importlib.resources from the package's equilibria/ data.
        # (In an editable install this maps back to the source equilibria/ dir;
        # in a built wheel it lives under site-packages/regcoil/equilibria/.)
        assert path.parent.name == "equilibria"


def test_aliases_return_identical_paths():
    a, b, c = examples("NCSX"), examples("li383"), examples("li383_1.4m")
    assert a.wout == b.wout == c.wout
    assert a.bnorm == b.bnorm == c.bnorm
    assert a.nescin == b.nescin == c.nescin


def test_dataset_is_frozen():
    ds = examples("NCSX")
    assert dataclasses.is_dataclass(ds)
    with pytest.raises(dataclasses.FrozenInstanceError):
        ds.wout = Path("/somewhere/else")  # type: ignore[misc]


@pytest.mark.parametrize("bad", ["nope", "", "   ", "li384", "ncsx-1.4m"])
def test_unknown_name_raises_keyerror(bad):
    with pytest.raises(KeyError):
        examples(bad)


def test_unknown_name_message_lists_available():
    with pytest.raises(KeyError) as excinfo:
        examples("does-not-exist")
    message = str(excinfo.value)
    assert "does-not-exist" in message
    for name in examples.available():
        assert name in message


def test_paths_feed_the_real_loaders():
    """The resolved paths are usable by the from_* loaders end to end."""
    ds = examples("NCSX")
    plasma = regcoil.PlasmaSurface.from_wout(ds.wout, ntheta=16, nzeta=16)
    plasma.set_bnormal_from_bnorm_file(ds.bnorm)
    coil = regcoil.CoilSurface.from_nescin(ds.nescin, nfp=plasma.nfp, ntheta=16, nzeta=16)
    assert plasma.nfp == 3
    assert coil.nfp == plasma.nfp
    assert len(coil.xm) > 0
