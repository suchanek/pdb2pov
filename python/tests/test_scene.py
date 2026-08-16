"""The scene writer: geometry, identifiers, and the shape of the output."""

from __future__ import annotations

import math

import pytest

from pypdb2pov import (
    Backdrop,
    Ground,
    ParseOptions,
    RadiiSet,
    SceneOptions,
    find_bonds,
    pov_identifier,
    prepare_structure,
    read_structure,
    scene_text,
)
from pypdb2pov.geometry import DEG_PER_RAD_CAM
from pypdb2pov.structure import Atom, Structure


def crambin(path: str, **scene_kwargs):
    structure, _ = read_structure(path, ParseOptions())
    options = SceneOptions(name="crambin", timestamp=False, **scene_kwargs)
    prepare_structure(structure, options)
    bonds = find_bonds(structure, options.bond_threshold) if options.ball_stick else None
    return structure, options, scene_text(structure, options, bonds)


# ----------------------------------------------------------------------
# Identifiers
# ----------------------------------------------------------------------


@pytest.mark.parametrize(
    "stem, expected",
    [
        ("crambin", "crambin"),
        ("out/1crn", "_1crn"),
        ("out/1crn.pov", "_1crn"),
        ("my molecule", "my_molecule"),
        ("a-b.c", "a_b_c"),
        ("1crn.inc", "_1crn"),
    ],
)
def test_identifiers_are_legal_povray(stem, expected):
    """`#declare out/1crn` does not parse; the C emitted it anyway."""
    assert pov_identifier(stem) == expected


# ----------------------------------------------------------------------
# Structure of the output
# ----------------------------------------------------------------------


def test_a_space_filling_scene_has_a_camera_and_no_bonds(crambin_pdb):
    _, _, text = crambin(crambin_pdb, backdrop=Backdrop.PLAIN)

    assert "camera {" in text
    assert "light_source {" in text
    assert "cylinder" not in text
    assert '#include "atoms_vdw.inc"' in text
    assert text.count("object { Atom_") == 327


def test_ball_and_stick_declares_a_bond_radius(crambin_pdb):
    _, _, text = crambin(crambin_pdb, ball_stick=True, bond_threshold=1.9)

    assert "#declare BOND_RADIUS = 0.17;" in text
    assert "#declare ATM_SCL = 0.30;" in text
    assert text.count("cylinder {") > 300


def test_object_only_writes_no_camera_and_restores_the_version(crambin_pdb):
    _, _, text = crambin(crambin_pdb, object_only=True)

    assert "camera {" not in text
    assert "light_source {" not in text
    assert text.startswith("//\n")
    assert "#declare crambin_pov_version = version;" in text
    assert text.rstrip().endswith("#version crambin_pov_version;")


def test_the_enclosing_radius_is_emitted_as_a_float(crambin_pdb):
    _, _, text = crambin(crambin_pdb)

    assert "#declare crambin_enclosing_radius = 18.759;" in text
    assert "//\tEnclosing Sphere: 18.759" in text


def test_the_camera_distance_quiltwright_documents_is_unchanged(crambin_pdb):
    """
    quiltwright's docs pin crambin at a camera distance of 40.075.  The
    truncated degree-to-radian constants that produce it are deliberate.
    """
    _, _, text = crambin(crambin_pdb)
    assert "  location  <0, 0, -40.075>" in text


def test_covalent_radii_select_the_other_include(crambin_pdb):
    _, _, text = crambin(crambin_pdb, radii=RadiiSet.COVALENT)
    assert '#include "atoms_covalent.inc"' in text
    assert '#include "atoms_vdw.inc"' not in text


def test_glass_atoms_merge_a_second_copy(crambin_pdb):
    _, _, text = crambin(crambin_pdb, ball_stick=True, glass_atoms=True, bond_threshold=1.9)

    assert "#declare crambin_obj_glass = merge {" in text
    assert "#declare ATM_SCL_B = 1.00;" in text
    assert "object { Atom_Glass_C" in text
    assert "  object { crambin_obj_glass }" in text


def test_backdrops_and_grounds(crambin_pdb):
    _, _, sky = crambin(crambin_pdb, backdrop=Backdrop.SKY, ground=Ground.CHECKER)
    assert "Gradient blue sky" in sky
    assert "checker color Black color White" in sky

    _, _, plain = crambin(crambin_pdb, backdrop=Backdrop.PLAIN, ground=Ground.PLAIN)
    assert "Plain white sky" in plain
    assert "color RichBlue" in plain


def test_unknown_elements_render_as_a_grey_sphere_by_default():
    """Before 2.1 an unrecognised element vanished from the scene."""
    structure = Structure(atoms=[Atom(0.0, 0.0, 0.0, name="U", element="U")])
    from pypdb2pov import assign_types

    assign_types(structure)

    options = SceneOptions(name="x", timestamp=False)
    prepare_structure(structure, options)
    text = scene_text(structure, options)

    assert "object { Atom_X" in text


def test_legacy_mode_drops_them_as_it_always_did():
    structure = Structure(atoms=[Atom(0.0, 0.0, 0.0, name="U", element="U")])
    from pypdb2pov import assign_types

    assign_types(structure, legacy=True)

    options = SceneOptions(name="x", timestamp=False, legacy_drop_unknown=True)
    prepare_structure(structure, options)
    text = scene_text(structure, options)

    assert "object { Atom_" not in text


def test_two_runs_of_the_same_conversion_are_byte_identical(crambin_pdb):
    _, _, first = crambin(crambin_pdb, ball_stick=True, bond_threshold=1.9)
    _, _, second = crambin(crambin_pdb, ball_stick=True, bond_threshold=1.9)
    assert first == second


# ----------------------------------------------------------------------
# Geometry
# ----------------------------------------------------------------------


def test_rotation_is_absolute_and_precedes_centring():
    structure = Structure(
        atoms=[Atom(1.0, 0.0, 0.0, element="C"), Atom(3.0, 0.0, 0.0, element="C")]
    )
    options = SceneOptions(zrot=90.0, name="x")
    prepare_structure(structure, options)

    # 90 degrees about z, then centred, then z flipped.  The tolerance is
    # 1e-3 rather than 1e-9 on purpose: the rotation uses the C's truncated
    # 57.29 for 180/pi, so "90 degrees" is 90.0009 degrees and always has
    # been.  Every existing scene was framed with that constant.
    # The pair spans x = 1..3, so after the rotation and centring the first
    # atom sits one angstrom off the centre along +y.
    assert structure.atoms[0].x == pytest.approx(0.0, abs=1e-3)
    assert structure.atoms[0].y == pytest.approx(1.0, abs=1e-3)
    assert structure.atoms[0].x != 0.0


def test_z_is_flipped_into_povrays_left_handed_world():
    structure = Structure(atoms=[Atom(0.0, 0.0, 1.0), Atom(0.0, 0.0, -1.0)])
    prepare_structure(structure, SceneOptions(name="x"))
    assert structure.atoms[0].z == -1.0


def test_the_camera_frames_the_widest_half_extent():
    structure = Structure(atoms=[Atom(-5.0, 0.0, 0.0), Atom(5.0, 0.0, 0.0)])
    options = SceneOptions(name="x", timestamp=False)
    prepare_structure(structure, options)
    text = scene_text(structure, options)

    extents = structure.extents()
    widest = max(
        (extents.xmax - extents.xmin) / 2.0,
        (extents.ymax - extents.ymin) / 2.0,
        (extents.zmax - extents.zmin) / 2.0,
    )
    expected = -widest / math.tan(22.0 / DEG_PER_RAD_CAM)
    assert f"location  <0, 0, {expected:.3f}>" in text
