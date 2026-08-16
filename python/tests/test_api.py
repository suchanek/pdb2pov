"""The library surface, which the C program never had."""

from __future__ import annotations

import os

import pytest

import pdb2pov
from pdb2pov import ParseOptions, SceneOptions, convert, read_structure


def test_convert_does_the_whole_pipeline(crambin_pdb, tmp_path):
    out = str(tmp_path / "crambin")
    structure, stats, path = convert(crambin_pdb, out, ball_stick=True, bond_threshold=1.9)

    assert path == out + ".pov"
    assert structure.natoms == 327
    assert stats.accepted == 327
    assert "cylinder {" in open(path).read()


def test_convert_takes_keywords_for_either_options_object(crambin_pdb, tmp_path):
    out = str(tmp_path / "chain_a")
    structure, _, _ = convert(crambin_pdb, out, chains="A", object_only=True)

    assert structure.natoms == 327  # crambin is one chain
    assert os.path.exists(out + ".inc")


def test_convert_rejects_an_unknown_keyword(crambin_pdb, tmp_path):
    with pytest.raises(TypeError, match="wobble"):
        convert(crambin_pdb, str(tmp_path / "x"), wobble=True)


def test_structure_is_a_sequence_of_atoms(crambin_pdb):
    structure, _ = read_structure(crambin_pdb, ParseOptions())

    assert len(structure) == 327
    assert structure[0].name == "N"
    assert [a.name for a in structure][:3] == ["N", "CA", "C"]


def test_structure_reports_what_it_holds(crambin_pdb):
    structure, _ = read_structure(crambin_pdb, ParseOptions())

    assert structure.chains() == ["A"]
    counts = structure.element_counts()
    assert counts["C"] == 202
    assert set(counts) == {"C", "N", "O", "S"}


def test_filtered_returns_a_new_structure(crambin_pdb):
    structure, _ = read_structure(crambin_pdb, ParseOptions())
    sulphurs = structure.filtered(lambda a: a.element == "S")

    assert sulphurs.natoms == 6
    assert structure.natoms == 327  # unchanged


def test_include_dir_ships_every_include_the_scenes_reference():
    directory = pdb2pov.include_dir()
    expected = {
        "atoms2.inc",
        "atoms_glass2.inc",
        "atoms_vdw.inc",
        "atoms_covalent.inc",
        "atoms_cpk.inc",
    }
    assert expected <= set(os.listdir(directory))


def test_every_element_in_the_table_has_its_declarations():
    """
    A row added to ELEMENTS without its POV-Ray declarations is a parse error
    at render time.  `make check` in the C tree catches that by rendering;
    this catches it without POV-Ray installed.
    """
    directory = pdb2pov.include_dir()
    solid = open(os.path.join(directory, "atoms2.inc")).read()
    glass = open(os.path.join(directory, "atoms_glass2.inc")).read()
    radii = {
        name: open(os.path.join(directory, name)).read()
        for name in ("atoms_vdw.inc", "atoms_covalent.inc", "atoms_cpk.inc")
    }

    for element in pdb2pov.ELEMENTS:
        assert f"#declare Atom_{element.pov_suffix} " in solid or (
            f"#declare Atom_{element.pov_suffix}\n" in solid
        ), element.symbol
        assert f"#declare Atom_Glass_{element.pov_suffix}" in glass, element.symbol
        for name, text in radii.items():
            assert f"#declare {element.symbol}_RAD" in text, f"{element.symbol} in {name}"


def test_the_public_names_are_all_importable():
    for name in pdb2pov.__all__:
        assert hasattr(pdb2pov, name), name


def test_scene_options_pick_the_right_include_file():
    from pdb2pov import RadiiSet

    assert RadiiSet.VDW.include_file == "atoms_vdw.inc"
    assert RadiiSet.COVALENT.include_file == "atoms_covalent.inc"
    assert SceneOptions().extension == ".pov"
    assert SceneOptions(object_only=True).extension == ".inc"
