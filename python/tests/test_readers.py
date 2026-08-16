"""Parsing: columns, filters, formats, and the things the C got wrong."""

from __future__ import annotations

import gzip
import os

import pytest

from pdb2pov import (
    ALL_MODELS,
    AltLocPolicy,
    InputFormat,
    ParseError,
    ParseOptions,
    read_pdb,
    read_structure,
)
from pdb2pov.readers import decode_hybrid36, detect_format, resolve_input_path


def read(path: str, **kwargs):
    return read_structure(path, ParseOptions(**kwargs))


# ----------------------------------------------------------------------
# Columns
# ----------------------------------------------------------------------


def test_reads_fixed_columns(data_dir):
    structure, stats = read(os.path.join(data_dir, "altloc.pdb"))

    first = structure.atoms[0]
    assert (first.name, first.res_name, first.chain_id) == ("CA", "ARG", "A")
    assert (first.x, first.y, first.z) == (0.0, 0.0, 0.0)
    assert first.element == "C"
    assert first.temp_factor == 10.0
    assert stats.freeform_fallback == 0


def test_free_form_fallback_is_accepted_and_reported(data_dir):
    """The C fell back to a free-form scan; dropping that would break files."""
    structure, stats = read(os.path.join(data_dir, "freeform.pdb"))

    assert structure.natoms == 2
    assert stats.freeform_fallback == 2
    assert structure.atoms[0].position == pytest.approx((13.5, 11.2, 8.4))


def test_short_records_do_not_read_past_the_end():
    structure, stats = read_pdb(["ATOM      1  N   CYS     3      13.507  11.238   8.398"])
    assert structure.natoms == 1
    assert stats.skipped_malformed == 0


def test_a_record_with_no_coordinates_is_reported_with_its_line(data_dir):
    lines = [
        "ATOM      1  N   CYS     3      13.507  11.238   8.398",
        "ATOM      2  CA  CYS     3      what   ever    nonsense",
    ]
    _, stats = read_pdb(lines)

    assert stats.skipped_malformed == 1
    assert stats.malformed_examples[0][0] == 2


def test_strict_mode_refuses_to_drop_a_record():
    lines = ["ATOM      2  CA  CYS     3      what   ever    nonsense"]
    with pytest.raises(ParseError):
        read_pdb(lines, ParseOptions(strict=True))


def test_long_lines_are_not_split(tmp_path):
    """The C read into a 256-byte buffer; a longer line became two records."""
    padded = (
        "ATOM      1  N   CYS     3      13.507  11.238   8.398  1.00  0.00"
        + " " * 400
        + "N"
    )
    structure, stats = read_pdb([padded])
    assert structure.natoms == 1
    assert stats.skipped_malformed == 0


# ----------------------------------------------------------------------
# Alternate conformations
# ----------------------------------------------------------------------


def test_altloc_default_keeps_blank_and_a(data_dir):
    structure, stats = read(os.path.join(data_dir, "altloc.pdb"))

    assert structure.natoms == 4
    assert stats.skipped_altloc == 3
    assert {a.alt_loc for a in structure.atoms} == {" ", "A"}


def test_altloc_all_keeps_every_conformer(data_dir):
    structure, _ = read(os.path.join(data_dir, "altloc.pdb"), altloc=AltLocPolicy.ALL)
    assert structure.natoms == 7


def test_altloc_first_matches_the_default_on_a_conventional_file(data_dir):
    default, _ = read(os.path.join(data_dir, "altloc.pdb"))
    first, _ = read(os.path.join(data_dir, "altloc.pdb"), altloc=AltLocPolicy.FIRST)
    assert [a.position for a in first] == [a.position for a in default]


def test_altloc_first_keeps_a_b_only_residue():
    """A residue labelled B and C with no A vanishes entirely under 'a'."""
    lines = [
        "ATOM      1  CA BARG A   1       1.500   0.000   0.000  0.50  0.00           C",
        "ATOM      2  CA CARG A   1       1.500   0.300   0.000  0.50  0.00           C",
    ]
    strict, _ = read_pdb(lines, ParseOptions(altloc=AltLocPolicy.A))
    lenient, _ = read_pdb(lines, ParseOptions(altloc=AltLocPolicy.FIRST))

    assert strict.natoms == 0
    assert lenient.natoms == 1
    assert lenient.atoms[0].alt_loc == "B"


def test_microheterogeneity_is_resolved_per_residue_not_per_atom(data_dir):
    """
    1CBN models residue 22 as both serine (altLocs A and B, 0.20 occupancy)
    and proline (altLoc C, 0.60), sharing one sequence position.  Choosing
    conformers atom by atom would take proline's ring and serine's hydroxyl
    from the same place and draw a residue that does not exist.
    """
    path = os.path.join(data_dir, "microhetero.pdb")

    for policy in AltLocPolicy:
        if policy is AltLocPolicy.ALL:
            continue
        structure, _ = read(path, altloc=policy)
        residues = {a.res_name for a in structure}
        assert len(residues) == 1, f"{policy.value} mixed {residues}"

    # 'a' and 'first' take the A conformer; 'occupancy' takes the major one.
    assert {a.res_name for a in read(path)[0]} == {"SER"}
    assert {a.res_name for a in read(path, altloc=AltLocPolicy.FIRST)[0]} == {"SER"}
    assert {a.res_name for a in read(path, altloc=AltLocPolicy.OCCUPANCY)[0]} == {"PRO"}


def test_altloc_occupancy_prefers_the_major_conformer(data_dir):
    path = os.path.join(data_dir, "occupancy.pdb")

    default, _ = read(path)
    major, _ = read(path, altloc=AltLocPolicy.OCCUPANCY)

    assert default.natoms == major.natoms == 2
    assert default.atoms[1].alt_loc == "A"  # the 0.30 conformer
    assert major.atoms[1].alt_loc == "B"  # the 0.70 one
    assert major.atoms[1].occupancy == 0.70


# ----------------------------------------------------------------------
# Models and filters
# ----------------------------------------------------------------------


def test_first_model_is_the_default(data_dir):
    structure, stats = read(os.path.join(data_dir, "models.pdb"))

    assert structure.natoms == 3
    assert {a.model for a in structure} == {1}
    assert stats.models_seen == [1, 2]
    assert stats.skipped_model == 3


def test_a_later_model_can_be_asked_for(data_dir):
    structure, _ = read(os.path.join(data_dir, "models.pdb"), model=2)
    assert structure.natoms == 3
    assert all(a.z == 2.0 for a in structure)


def test_all_models_can_be_kept(data_dir):
    structure, _ = read(os.path.join(data_dir, "models.pdb"), model=ALL_MODELS)
    assert structure.natoms == 6


def test_chain_filter(data_dir):
    structure, stats = read(os.path.join(data_dir, "models.pdb"), chains="B")
    assert structure.natoms == 1
    assert stats.skipped_chain == 2


def test_hetatm_and_water_filters(data_dir):
    path = os.path.join(data_dir, "haem.pdb")

    everything, _ = read(path)
    protein_only, stats = read(path, keep_hetatm=False)

    assert everything.natoms == 5
    assert protein_only.natoms == 1
    assert stats.skipped_hetatm == 4


# ----------------------------------------------------------------------
# Elements
# ----------------------------------------------------------------------


def test_element_column_wins_over_the_name(data_dir):
    structure, stats = read(os.path.join(data_dir, "altloc.pdb"))
    assert [a.element for a in structure] == ["C"] * 4
    assert stats.no_element_column == 0


def test_inference_reads_haem_nitrogens_as_nitrogen(data_dir):
    """`NA` in a haem is N-alpha.  The pre-2.1 guess called it sodium."""
    structure, stats = read(os.path.join(data_dir, "haem.pdb"))
    by_name = {a.name: a.element for a in structure}

    assert by_name["NA"] == "N"
    assert by_name["NB"] == "N"
    assert by_name["FE"] == "FE"  # right-justified into columns 13-14
    assert by_name["CD"] == "C"  # glutamate C-delta, not cadmium
    assert by_name["ZN"] == "ZN"  # atom name agrees with the residue name
    assert stats.no_element_column == 5
    assert stats.inferred_elements == 5


def test_legacy_elements_reproduce_the_old_mistakes(data_dir):
    from pdb2pov import element_symbol

    structure, _ = read(os.path.join(data_dir, "haem.pdb"), legacy_elements=True)
    by_name = {a.name: element_symbol(a.type_index) for a in structure}

    assert by_name["NA"] == "N"
    assert by_name["FE"] == "FE"
    assert by_name["ZN"] == "X"  # the guess knew seven letters; Z was not one


def test_ambiguous_names_are_reported_not_silent():
    """A hand-edited ' FE ' reads as fluorine per the format, and says so."""
    lines = ["HETATM   99  FE  XXX     4      33.265  10.443  10.307"]
    structure, stats = read_pdb(lines)

    assert structure.atoms[0].element == "F"
    assert stats.ambiguous_names == [("FE", "F", "FE")]
    assert any("two-letter element" in line for line in stats.lines())


def test_deuterium_draws_as_hydrogen():
    from pdb2pov import ELEMENT_UNKNOWN, element_symbol

    lines = [
        "ATOM      1  D1  GLY A   1       0.000   0.000   0.000  1.00  0.00           D"
    ]
    structure, _ = read_pdb(lines)
    from pdb2pov import assign_types

    assign_types(structure)
    assert structure.atoms[0].type_index != ELEMENT_UNKNOWN
    assert element_symbol(structure.atoms[0].type_index) == "H"


# ----------------------------------------------------------------------
# Formats and paths
# ----------------------------------------------------------------------


def test_hybrid36_serial_numbers():
    assert decode_hybrid36("99999", 5) == 99999
    assert decode_hybrid36("A0000", 5) == 100000
    assert decode_hybrid36("A0001", 5) == 100001
    assert decode_hybrid36("a0000", 5) == 43770016
    assert decode_hybrid36("   12", 5) == 12
    assert decode_hybrid36("     ", 5) == 0


def test_a_big_serial_does_not_make_a_record_look_broken():
    line = "ATOM  A0000  N   CYS     3      13.507  11.238   8.398"
    structure, stats = read_pdb([line])
    assert stats.skipped_malformed == 0
    assert structure.atoms[0].serial == 100000


def test_gzip_is_read_transparently(tmp_path, data_dir):
    source = open(os.path.join(data_dir, "altloc.pdb"), "rb").read()
    packed = tmp_path / "altloc.pdb.gz"
    with gzip.open(packed, "wb") as handle:
        handle.write(source)

    structure, _ = read_structure(str(packed))
    assert structure.natoms == 4


def test_a_bare_stem_still_finds_the_pdb(data_dir):
    stem = os.path.join(data_dir, "altloc")
    assert resolve_input_path(stem).endswith("altloc.pdb")

    structure, _ = read_structure(stem)
    assert structure.natoms == 4


def test_format_detection_prefers_content_over_extension(tmp_path, data_dir):
    misnamed = tmp_path / "actually.pdb"
    misnamed.write_text(open(os.path.join(data_dir, "mini.cif")).read())

    assert detect_format(str(misnamed), open(misnamed).readlines()) is InputFormat.CIF

    structure, _ = read_structure(str(misnamed))
    assert structure.source_format == "cif"
    assert structure.natoms == 4
