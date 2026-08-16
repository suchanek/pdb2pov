"""The mmCIF reader, which is the format the C could never read at all."""

from __future__ import annotations

import os

from pypdb2pov import ALL_MODELS, ParseOptions, read_structure
from pypdb2pov.mmcif import read_mmcif, tokenise


def test_reads_the_atom_site_loop(data_dir):
    structure, stats = read_structure(os.path.join(data_dir, "mini.cif"))

    assert structure.source_format == "cif"
    assert structure.natoms == 4  # model 1 only
    assert [a.element for a in structure] == ["N", "C", "O", "ZN"]
    assert [a.name for a in structure] == ["N", "CA", "O5'", "ZN"]
    assert stats.models_seen == [1, 2]


def test_quoted_values_survive_the_tokeniser(data_dir):
    structure, _ = read_structure(os.path.join(data_dir, "mini.cif"))
    assert structure.atoms[2].name == "O5'"


def test_semicolon_text_blocks_are_one_token(data_dir):
    structure, _ = read_structure(os.path.join(data_dir, "mini.cif"))
    assert structure.title.startswith("A four-atom fragment")
    assert "\n" in structure.title


def test_auth_columns_are_preferred(data_dir):
    """--chain A should mean the chain the literature calls A."""
    structure, _ = read_structure(
        os.path.join(data_dir, "mini.cif"), ParseOptions(chains="B")
    )
    assert structure.natoms == 1
    assert structure.atoms[0].element == "ZN"


def test_models(data_dir):
    path = os.path.join(data_dir, "mini.cif")

    second, _ = read_structure(path, ParseOptions(model=2))
    every, _ = read_structure(path, ParseOptions(model=ALL_MODELS))

    assert second.natoms == 2
    assert all(a.z == 9.0 for a in second)
    assert every.natoms == 6


def test_formal_charge_is_read(data_dir):
    structure, _ = read_structure(os.path.join(data_dir, "mini.cif"))
    assert structure.atoms[3].charge == 2.0


def test_a_loop_that_is_not_atom_site_is_skipped():
    text = """\
data_X
loop_
_chem_comp.id
_chem_comp.type
GLY 'peptide linking'
ALA 'peptide linking'
loop_
_atom_site.group_PDB
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
ATOM C CA 1.0 2.0 3.0
"""
    structure, _ = read_mmcif(text.splitlines())
    assert structure.natoms == 1
    assert structure.atoms[0].position == (1.0, 2.0, 3.0)


def test_null_values_become_empty():
    tokens = list(tokenise(["ATOM . ? 'a b' 1.0"]))
    assert [t[0] for t in tokens] == ["ATOM", ".", "?", "a b", "1.0"]
    assert [t[1] for t in tokens] == [False, False, False, True, False]


def test_a_key_with_no_value_does_not_swallow_the_next_keyword():
    """
    A truncated ``_struct.title`` immediately followed by ``loop_`` used to
    consume the keyword, hiding every atom behind it.
    """
    text = """\
data_X
_struct.title
loop_
_atom_site.type_symbol
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
C 1.0 2.0 3.0
"""
    structure, _ = read_mmcif(text.splitlines())
    assert structure.natoms == 1
    assert structure.title == "X" or structure.title == ""


# ----------------------------------------------------------------------
# Against a real wwPDB file
# ----------------------------------------------------------------------


def test_the_real_1crn_atom_site_loop_matches_the_pdb_file(data_dir, crambin_pdb):
    """
    The fixture is the ``_atom_site`` loop of the wwPDB's own ``1CRN.cif``.
    Reading it must give exactly what reading ``1CRN.pdb`` gives -- same
    atoms, same order, same numbers -- because they are the same entry.
    """
    from pypdb2pov import ParseOptions

    pdb, _ = read_structure(crambin_pdb, ParseOptions())
    cif, _ = read_structure(os.path.join(data_dir, "1crn_atom_site.cif"))

    assert cif.natoms == pdb.natoms == 327
    assert [a.name for a in cif] == [a.name for a in pdb]
    assert [a.element for a in cif] == [a.element for a in pdb]
    assert [a.position for a in cif] == [a.position for a in pdb]
    assert [a.res_seq for a in cif] == [a.res_seq for a in pdb]
    assert [a.temp_factor for a in cif] == [a.temp_factor for a in pdb]
    assert cif.title.startswith("WATER STRUCTURE OF A HYDROPHOBIC PROTEIN")


def test_the_same_entry_gives_the_same_scene_either_way(data_dir, crambin_pdb):
    """PDB in or mmCIF in, the scene that comes out is the same file."""
    from pypdb2pov import (
        ParseOptions,
        SceneOptions,
        find_bonds,
        prepare_structure,
        scene_text,
    )

    def render(path: str) -> str:
        structure, _ = read_structure(path, ParseOptions())
        options = SceneOptions(
            name="crambin", timestamp=False, ball_stick=True, bond_threshold=1.9
        )
        prepare_structure(structure, options)
        return scene_text(structure, options, find_bonds(structure, 1.9))

    from_pdb = render(crambin_pdb)
    from_cif = render(os.path.join(data_dir, "1crn_atom_site.cif"))

    # Only the source file named on the "Prepared by" line may differ.
    strip = lambda text: [ln for ln in text.splitlines() if "Prepared by" not in ln]
    assert strip(from_pdb) == strip(from_cif)


def test_a_truncated_final_row_is_reported_not_silently_kept():
    text = """\
data_X
loop_
_atom_site.type_symbol
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
C 1.0 2.0 3.0
N 4.0
"""
    structure, stats = read_mmcif(text.splitlines())
    assert structure.natoms == 1
    assert stats.skipped_malformed == 1
