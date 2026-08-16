"""Shared fixtures: the test data directory and the repository it lives in."""

from __future__ import annotations

import os
import subprocess

import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "data")

#: The C sources sit one directory up from the Python package.  The tests
#: that compare the two are skipped when they are not there, so the package
#: still tests clean after being copied somewhere else -- which is the point
#: of it being a stand-alone directory.
REPO_ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))


@pytest.fixture(scope="session")
def data_dir() -> str:
    return DATA


@pytest.fixture(scope="session")
def repo_root() -> str:
    return REPO_ROOT


@pytest.fixture(scope="session")
def crambin_pdb() -> str:
    path = os.path.join(REPO_ROOT, "1CRN.pdb")
    if not os.path.exists(path):
        pytest.skip("1CRN.pdb is not next to the package")
    return path


@pytest.fixture(scope="session")
def c_pdb2pov(tmp_path_factory) -> str:
    """
    Build the C program and return its path, or skip.

    Built into a temporary directory so the tests never disturb whatever the
    developer has in the source tree.
    """
    source = os.path.join(REPO_ROOT, "pdb2pov.c")
    if not os.path.exists(source):
        pytest.skip("the C sources are not next to the package")

    binary = os.path.join(REPO_ROOT, "pdb2pov")
    if os.path.exists(binary) and os.access(binary, os.X_OK):
        return binary

    out = str(tmp_path_factory.mktemp("cbuild") / "pdb2pov")
    try:
        subprocess.run(
            ["cc", "-std=c17", "-O1", "-o", out, source,
             os.path.join(REPO_ROOT, "util.c"), "-I", REPO_ROOT, "-lm"],
            check=True,
            capture_output=True,
        )
    except (OSError, subprocess.CalledProcessError) as exc:  # pragma: no cover
        pytest.skip(f"cannot build the C reference: {exc}")
    return out
