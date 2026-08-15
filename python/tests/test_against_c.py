"""
Differential tests against ``pdb2pov.c``.

The port's promise is that an existing command line produces the same scene,
byte for byte, apart from the header's "Prepared by" line -- which carries the
timestamp and so differs between any two runs of either program.  These tests
build the C program and diff the two.

They skip when the C sources are not alongside the package, so the Python
directory still tests clean after being copied elsewhere.
"""

from __future__ import annotations

import os
import subprocess
import sys

import pytest

FLAG_SETS = [
    ["-v", "-p"],
    ["-c", "-b", "-d", "1.9"],
    ["-b", "-d", "1.9", "-o"],
    ["-q", "-d", "1.9", "-s", "-h"],
    ["-b", "-d", "2.2", "-g", "-a"],
    ["-v", "-s", "-x", "90", "-y", "30", "-z", "12", "-r", "1.4"],
    ["-b", "-d", "1.9", "--chain", "A"],
    ["-b", "-d", "1.9", "--legacy-elements"],
    ["-b", "-d", "1.9", "--keep-altlocs"],
]


def normalise(path: str, stem: str) -> list[str]:
    """
    Strip what is allowed to differ: the timestamp line, and the output stem,
    which appears in every declared identifier and differs because the two
    programs write to different files.
    """
    with open(path, encoding="utf-8") as handle:
        return [
            line.replace(stem, "MOLECULE")
            for line in handle
            if "Prepared by" not in line
        ]


def run_c(binary: str, source: str, out_stem: str, flags: list[str]) -> int:
    return subprocess.run(
        [binary, source, out_stem, *flags], capture_output=True, text=True
    ).returncode


def run_python(source: str, out_stem: str, flags: list[str]) -> int:
    env = dict(os.environ)
    root = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src"))
    env["PYTHONPATH"] = root + os.pathsep + env.get("PYTHONPATH", "")
    return subprocess.run(
        [sys.executable, "-m", "pdb2pov", source, out_stem, "--quiet", *flags],
        capture_output=True,
        text=True,
        env=env,
    ).returncode


def produced(stem: str) -> str | None:
    for extension in (".pov", ".inc"):
        if os.path.exists(stem + extension):
            return stem + extension
    return None


@pytest.mark.parametrize("flags", FLAG_SETS, ids=lambda f: " ".join(f))
def test_crambin_scene_matches_the_c_program(c_pdb2pov, crambin_pdb, tmp_path, flags):
    source = crambin_pdb[: -len(".pdb")]
    c_stem = str(tmp_path / "from_c")
    py_stem = str(tmp_path / "from_py")

    assert run_c(c_pdb2pov, source, c_stem, flags) == 0
    assert run_python(source, py_stem, flags) == 0

    c_out, py_out = produced(c_stem), produced(py_stem)
    assert c_out and py_out

    assert normalise(c_out, "from_c") == normalise(py_out, "from_py")


def test_the_element_grid_matches(c_pdb2pov, repo_root, tmp_path):
    """
    One atom of every element in the table.  This is the file `make check`
    renders to prove the element table and the include files agree; the port
    has to place all thirty-three the same way.
    """
    source = os.path.join(repo_root, "elements")
    if not os.path.exists(source + ".pdb"):
        pytest.skip("elements.pdb is not next to the package")

    c_stem = str(tmp_path / "from_c")
    py_stem = str(tmp_path / "from_py")

    assert run_c(c_pdb2pov, source, c_stem, ["-v", "-p"]) == 0
    assert run_python(source, py_stem, ["-v", "-p"]) == 0

    assert normalise(c_stem + ".pov", "from_c") == normalise(py_stem + ".pov", "from_py")

    # Atom_X is the generic grey sphere; Atom_Xe is xenon, so the boundary
    # matters.  Nothing in the table may fall through to the generic.
    import re

    assert not re.search(r"\bAtom_X\b", open(py_stem + ".pov").read())


def test_legacy_mode_reproduces_the_old_name_guessing(c_pdb2pov, repo_root, tmp_path):
    """
    ``test.pdb`` has no element column, so it is the file where the improved
    inference and the 1993 guess disagree.  Under ``--legacy-elements`` they
    must not.
    """
    source = os.path.join(repo_root, "test")
    if not os.path.exists(source + ".pdb"):
        pytest.skip("test.pdb is not next to the package")

    c_stem = str(tmp_path / "from_c")
    py_stem = str(tmp_path / "from_py")
    flags = ["-b", "-d", "1.9", "-p", "--legacy-elements"]

    assert run_c(c_pdb2pov, source, c_stem, flags) == 0
    assert run_python(source, py_stem, flags) == 0

    assert normalise(c_stem + ".pov", "from_c") == normalise(py_stem + ".pov", "from_py")


def test_the_bundled_includes_are_the_same_files(repo_root):
    """
    The package ships its own copy of the POV-Ray includes so it can be used
    on its own.  Copies drift; this is what notices.
    """
    from pdb2pov import include_dir

    for name in os.listdir(include_dir()):
        upstream = os.path.join(repo_root, name)
        if not os.path.exists(upstream):
            pytest.skip(f"{name} is not next to the package")
        with open(os.path.join(include_dir(), name), "rb") as ours, open(
            upstream, "rb"
        ) as theirs:
            assert ours.read() == theirs.read(), f"{name} has drifted from the C tree's copy"


def test_the_element_tables_agree(repo_root):
    """The C's ELEMENTS array and the Python one must list the same rows."""
    source = os.path.join(repo_root, "pdb2pov.c")
    if not os.path.exists(source):
        pytest.skip("the C sources are not next to the package")

    import re

    from pdb2pov import ELEMENTS

    text = open(source, encoding="utf-8").read()
    body = text.split("static const Element ELEMENTS[] = {", 1)[1].split("};", 1)[0]
    rows = re.findall(r'\{\s*"([A-Za-z]+)"\s*,\s*"([A-Za-z]+)"\s*\}', body)

    assert rows == [(e.symbol, e.pov_suffix) for e in ELEMENTS]


def test_the_radii_match_the_include_files(repo_root):
    """
    The covalent bonding mode reads radii from a table in Python rather than
    from the POV-Ray includes.  They have to be the same numbers.
    """
    import re

    from pdb2pov.elements import COVALENT_RADII, VDW_RADII

    for filename, table in (
        ("atoms_vdw.inc", VDW_RADII),
        ("atoms_covalent.inc", COVALENT_RADII),
    ):
        path = os.path.join(repo_root, filename)
        if not os.path.exists(path):
            pytest.skip(f"{filename} is not next to the package")

        declared = {
            symbol: float(value)
            for symbol, value in re.findall(
                r"#declare\s+([A-Z]+)_RAD\s*=\s*([0-9.]+)\s*;", open(path).read()
            )
            if symbol != "X"
        }
        assert declared == table, f"{filename} disagrees with the Python table"
