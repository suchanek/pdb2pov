"""
The POV-Ray scene writer.

Output is byte-identical to ``pdb2pov.c`` 2.2 for the same input and flags,
apart from the timestamp in the header comment.  That is a deliberate
constraint rather than a coincidence: existing scenes, and quiltwright's
documented camera distances and enclosing radii, are derived from these exact
numbers, so the port has to reproduce the arithmetic and the formatting
rather than tidy either up.  ``tests/test_against_c.py`` diffs the two
programs' output to keep it that way.

Scenes are POV-Ray 3.7 and depend on the bundled include files
(``atoms2.inc`` and the radius sets) being on POV-Ray's library path.  See
:func:`pypdb2pov.include_dir`.
"""

from __future__ import annotations

import math
import os
import time
from typing import IO, Sequence

from .elements import ELEMENT_UNKNOWN, element_pov_suffix
from .geometry import DEG_PER_RAD_CAM, rotate_structure
from .options import Backdrop, Ground, SceneOptions
from .structure import Extents, Structure

__all__ = [
    "PYPDB2POV_VERSION",
    "PDB2POV_VERSION",
    "POV_VERSION",
    "BOND_RAD",
    "SPHERE_FUDGE",
    "pov_identifier",
    "prepare_structure",
    "scene_text",
    "write_scene",
]

#: This package's own version.  It moves independently of the number below,
#: so a port-only fix does not have to claim the C program changed.
PYPDB2POV_VERSION = "0.1.0"

#: The pdb2pov release these scenes are those of, kept in step with
#: ``PDB2POV_VERSION`` in ``pdb2pov.h``.  The scenes are byte-identical to
#: that release's, and a header claiming otherwise would be misleading, so
#: this tracks the C rather than the package.
PDB2POV_VERSION = "2.2"

#: POV-Ray language level the emitted scenes are written against.
POV_VERSION = "3.7"

#: Radius of the cylinders drawn between bonded atoms, in angstroms.
BOND_RAD = 0.17

#: The bounding sphere is grown by this fraction so it clears the outermost atom.
SPHERE_FUDGE = 0.02

CAMERA_HALF_ANGLE_DEG = 22.0
LIGHT_HALF_ANGLE_DEG = 20.0

_TOOL = "pypdb2pov"


def pov_identifier(stem: str, fallback: str = "molecule") -> str:
    """
    Derive a legal POV-Ray identifier from an output path.

    POV identifiers are alphanumerics and underscores and may not begin with
    a digit, so a stem like ``out/1crn`` -- perfectly good as a path -- needs
    both the directory dropped and a leading underscore added.  A trailing
    ``.pov`` or ``.inc`` is dropped too, which the C never had to do because
    it only ever took stems.
    """
    base = os.path.basename(stem)
    for suffix in (".pov", ".inc"):
        if base.lower().endswith(suffix):
            base = base[: -len(suffix)]
            break

    if not base:
        base = fallback
    if base[0].isdigit():
        base = "_" + base

    return "".join(ch if (ch.isalnum() and ch.isascii()) or ch == "_" else "_" for ch in base)


def prepare_structure(structure: Structure, options: SceneOptions) -> tuple[float, float, float]:
    """
    Rotate, centre and flip a structure into POV-Ray's world, in place.

    Returns the centring shift.  This is the C's ``main()`` sequence, and the
    order matters: rotation is about the original origin, not the centroid.
    """
    rotate_structure(structure, options.xrot, options.yrot, options.zrot)
    shift = structure.center()
    structure.flip_z()
    return shift


def scene_text(
    structure: Structure,
    options: SceneOptions,
    bonds: Sequence[tuple[int, int]] | None = None,
) -> str:
    """
    Render the scene as a string.

    ``structure`` must already have been through :func:`prepare_structure`
    and had its element types assigned.  ``bonds`` is ignored unless
    :attr:`SceneOptions.ball_stick` is set.
    """
    out: list[str] = []
    extents = structure.extents()
    sphere_rad = structure.enclosing_radius()
    sphere_rad += SPHERE_FUDGE * sphere_rad

    atom_scale = options.radii_scale
    if options.ball_stick:
        atom_scale *= 0.3

    ident = options.name

    _write_header(out, structure, extents, sphere_rad, options)

    if options.object_only:
        # An include has to leave the language version as it found it, so
        # the host scene is not silently switched to 3.7 from whatever it set.
        out.append(f"#declare {ident}_pov_version = version;\n")
        out.append(f"#version {POV_VERSION};\n\n")
    else:
        out.append(f"#version {POV_VERSION};\n\n")
        # Required from 3.7 onwards; POV-Ray errors out without it.
        out.append("global_settings { assumed_gamma 1.0 }\n\n")

        out.append('#include "colors.inc"\n')
        out.append('#include "shapes.inc"\n')
        out.append('#include "textures.inc"\n')
        out.append(f'#include "{options.radii.include_file}"\n')
        out.append('#include "atoms2.inc"\n')
        if options.glass_atoms:
            out.append('#include "atoms_glass2.inc"\n')
        out.append("\n")

        _write_camera(out, extents)
        _write_light(out, extents, options.area_light)

        if options.backdrop is Backdrop.SKY:
            _write_sky(out)
        elif options.backdrop is Backdrop.PLAIN:
            _write_sky_plain(out)

        if options.ground is Ground.PLAIN:
            _write_ground(out, extents)
        elif options.ground is Ground.CHECKER:
            _write_check(out, extents)

    # Float declarations need a terminating semicolon from POV-Ray 3.5
    # onwards.  Object and texture declarations do not, and do not get one.
    if options.ball_stick:
        out.append(f"#declare BOND_RADIUS = {BOND_RAD:.2f};\n")
    if options.glass_atoms:
        out.append(f"#declare ATM_SCL_B = {atom_scale / 0.3:.2f};\n")
    out.append(f"#declare ATM_SCL = {atom_scale:.2f};\n\n")

    # The glass shell, if asked for, is merged so interior faces vanish.
    if options.glass_atoms:
        out.append(f"#declare {ident}_obj_glass = merge {{\n")
        _write_atoms(out, structure, options, "ATM_SCL_B", glass=True)
        out.append("}\n\n")

    out.append(f"#declare {ident}_obj = union {{\n")
    _write_atoms(out, structure, options, "ATM_SCL", glass=False)

    if options.ball_stick and bonds:
        atoms = structure.atoms
        for a, b in bonds:
            first, second = atoms[a], atoms[b]
            out.append(
                f"  cylinder {{ <{first.x:.3f}, {first.y:.3f}, {first.z:.3f}>, "
                f"<{second.x:.3f}, {second.y:.3f}, {second.z:.3f}>, "
                "BOND_RADIUS texture { White_Bond } }\n"
            )

    if options.glass_atoms:
        out.append(f"  object {{ {ident}_obj_glass }}\n")

    out.append("}\n\n")

    # The enclosing sphere radius as a float a host scene can read.  The 1993
    # code spent it on bounded_by { sphere { ... } }, which POV-Ray 3.x
    # rejects as redundant; the number is still what sets the depth budget
    # when these scenes are rendered for a holographic display.
    out.append(f"#declare {ident}_enclosing_radius = {sphere_rad:.3f};\n")
    out.append(f"#declare {ident} = object {{ {ident}_obj }}\n")

    if options.object_only:
        out.append(f"\n#version {ident}_pov_version;\n")
    else:
        out.append(f"\nobject {{ {ident} }}\n")

    return "".join(out)


def write_scene(
    structure: Structure,
    options: SceneOptions,
    destination: str | IO[str],
    bonds: Sequence[tuple[int, int]] | None = None,
) -> str:
    """
    Write the scene to a path or an open stream, returning what was written.
    """
    text = scene_text(structure, options, bonds)
    if isinstance(destination, str):
        with open(destination, "w", encoding="ascii", newline="\n") as handle:
            handle.write(text)
    else:
        destination.write(text)
    return text


# ----------------------------------------------------------------------
# Sections
# ----------------------------------------------------------------------


def _write_header(
    out: list[str],
    structure: Structure,
    e: Extents,
    sphere_rad: float,
    options: SceneOptions,
) -> None:
    stamp = time.strftime("%Y-%m-%d %H:%M:%S") if options.timestamp else "an unrecorded date"
    origin = f" from {structure.source}" if structure.source else ""

    # Everything that varies between two runs of the same conversion lives on
    # this one line, so the rest of the file can be compared byte for byte --
    # both against a previous run and against the C program.
    out.append("//\n")
    out.append(
        f"// Prepared by {_TOOL} {PYPDB2POV_VERSION} "
        f"(pdb2pov {PDB2POV_VERSION}){origin} on {stamp}\n"
    )
    out.append("// Author: Eric G. Suchanek, Ph.D.\n")
    out.append("//\n")
    out.append(f"//\tAtoms: {structure.natoms:4d}\n")
    out.append(
        f"//\tExtent:\tXmin: {e.xmin:.3f} Xmax: {e.xmax:.3f},\n"
        f"//\t\tYmin: {e.ymin:.3f}, Ymax: {e.ymax:.3f}\n"
    )
    out.append(f"//\t\tZmin: {e.zmin:.3f} Zmax: {e.zmax:.3f}\n")
    out.append(f"//\tEnclosing Sphere: {sphere_rad:.3f}\n")
    out.append("//\n\n")


def _write_camera(out: list[str], e: Extents) -> None:
    """
    Camera placement.

    The molecule is centred on the origin, so the camera sits back along -z
    far enough that the widest half-extent subtends
    :data:`CAMERA_HALF_ANGLE_DEG`.  ``right x*image_width/image_height``
    replaces the 2.x fixed ``right <4/3,0,0>``: at 4:3 the two are identical,
    at any other aspect the fixed form rendered non-square pixels.  The
    vertical field of view is unchanged at 2*atan(0.5) = 53.13 degrees.
    """
    xavg = (e.xmax - e.xmin) / 2.0
    yavg = (e.ymax - e.ymin) / 2.0
    zavg = (e.zmax - e.zmin) / 2.0
    theta = CAMERA_HALF_ANGLE_DEG / DEG_PER_RAD_CAM

    dist = min(-xavg / math.tan(theta), -yavg / math.tan(theta), -zavg / math.tan(theta))

    out.append("camera {\n")
    out.append(f"  location  <0, 0, {dist:.3f}>\n")
    out.append("  look_at   <0, 0, 0>\n")
    out.append("  direction <0, 0, 1>\n")
    out.append("  up        <0, 1, 0>\n")
    out.append("  right     x*image_width/image_height\n")
    out.append("}\n\n")


def _write_light(out: list[str], e: Extents, area_light: bool) -> None:
    """
    A single light off the upper corner.  In POV-Ray 3.x a light_source is
    not an object, so the 2.x ``object { light_source { ... } }`` is gone.
    """
    xavg = (e.xmax - e.xmin) / 2.0
    theta = LIGHT_HALF_ANGLE_DEG / DEG_PER_RAD_CAM
    dist = -xavg / math.tan(theta)

    out.append("light_source {\n")
    out.append(f"  <{e.xmax:.3f}, {e.ymax:.3f}, {dist:.3f}>\n")
    out.append("  color White\n")
    if area_light:
        out.append(f"  area_light <{e.xmax / 2.0:.3f}, 0, 0>, <0, 0, {e.xmax / 2.0:.3f}>, 9, 9\n")
        out.append("  adaptive 1\n")
        out.append("  jitter\n")
    out.append("}\n\n")


def _write_sky(out: list[str]) -> None:
    out.append("// Gradient blue sky with white clouds\n")
    out.append("sphere {\n")
    out.append("  <0, 0, 0>, 3000\n")
    out.append("  pigment {\n")
    out.append("    gradient y\n")
    out.append("    color_map {\n")
    out.append("      [0.0 rgb <0.30, 0.30, 1.00>]\n")
    out.append("      [0.8 rgb <0.70, 0.70, 1.00>]\n")
    out.append("      [1.0 rgb <0.90, 0.90, 1.00>]\n")
    out.append("    }\n")
    out.append("    scale <3000, 3000, 3000>\n")
    out.append("    quick_color rgb <0.70, 0.70, 1.00>\n")
    out.append("  }\n")
    out.append("  finish {\n")
    out.append("    ambient 0.7\n")
    out.append("    diffuse 0  // clouds must not shadow the sky\n")
    out.append("  }\n")
    out.append("}\n\n")

    out.append("sphere {\n")
    out.append("  <0, 0, 0>, 2590\n")
    out.append("  pigment {\n")
    out.append("    bozo\n")
    out.append("    turbulence 0.5\n")
    out.append("    color_map {\n")
    out.append("      [0.0 rgbf <1.00, 1.00, 1.00, 1.00>]\n")
    out.append("      [0.6 rgbf <1.00, 1.00, 1.00, 1.00>]\n")
    out.append("      [0.8 rgb  <1.00, 1.00, 1.00>]\n")
    out.append("      [1.0 rgb  <0.80, 0.80, 0.80>]\n")
    out.append("    }\n")
    out.append("    quick_color rgb <0.70, 0.70, 1.00>\n")
    out.append("    scale <1000, 200, 1000>\n")
    out.append("  }\n")
    out.append("  finish { ambient 0.7 diffuse 0 }\n")
    out.append("}\n\n")


def _write_sky_plain(out: list[str]) -> None:
    out.append("// Plain white sky\n")
    out.append("sphere {\n")
    out.append("  <0, 0, 0>, 3000\n")
    out.append("  pigment { color White }\n")
    out.append("  finish { ambient 0.7 diffuse 0 }\n")
    out.append("}\n\n")


def _ground_level(e: Extents) -> float:
    """Y of the ground plane, set below the lowest atom."""
    yavg = (e.ymax - e.ymin) / 2.0
    return -1.0 * (abs(e.ymin) + abs(yavg) / 2.0)


def _write_ground(out: list[str], e: Extents) -> None:
    out.append("plane {\n")
    out.append(f"  y, {_ground_level(e):.3f}\n")
    out.append("  pigment { color RichBlue quick_color RichBlue }\n")
    out.append("  normal { bumps 0.2 }\n")
    out.append("}\n\n")


def _write_check(out: list[str], e: Extents) -> None:
    xavg = (e.xmax - e.xmin) / 2.0
    yavg = (e.ymax - e.ymin) / 2.0

    out.append("// The beloved raytracer checkerboard\n")
    out.append("plane {\n")
    out.append(f"  y, {_ground_level(e):.3f}\n")
    out.append("  pigment {\n")
    out.append("    checker color Black color White\n")
    out.append(f"    scale {(xavg + yavg) / 3.0:.3f}\n")
    out.append("  }\n")
    out.append("  finish { ambient 0.2 diffuse 0.8 }\n")
    out.append("}\n\n")


def _write_atoms(
    out: list[str],
    structure: Structure,
    options: SceneOptions,
    scale_id: str,
    glass: bool,
) -> None:
    """
    One merge/union member per atom.

    Under legacy element handling an atom with no dedicated texture is
    dropped, which is what happened before 2.1; reproducing an old render
    means reproducing the omission.  By default such atoms render as the
    generic ``Atom_X`` and nothing disappears.
    """
    prefix = "Atom_Glass_" if glass else "Atom_"
    legacy_drop = options.legacy_drop_unknown

    for atom in structure.atoms:
        if legacy_drop and atom.type_index == ELEMENT_UNKNOWN:
            continue
        name = prefix + element_pov_suffix(atom.type_index)
        out.append(
            f"  object {{ {name} scale <{scale_id}, {scale_id}, {scale_id}> "
            f"translate <{atom.x:.3f}, {atom.y:.3f}, {atom.z:.3f}> }}\n"
        )
