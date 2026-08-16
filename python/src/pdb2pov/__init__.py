"""
pdb2pov -- convert Brookhaven PDB and PDBx/mmCIF structures into POV-Ray
scenes.

A Python port of ``pdb2pov.c`` 2.2 by Eric G. Suchanek, Ph.D., which has been
turning atomic coordinates into ray-traced pictures since 1993.  The scenes
this writes are byte-identical to the C program's for the same input and
flags; what is new is everything around the scene writer -- an mmCIF reader,
compressed input, model and chain selection, element inference that knows
more than seven first letters, a linear-time bond search, and a library API
rather than only a command line.

Typical use::

    from pdb2pov import convert

    convert("1crn.pdb", "crambin", ball_stick=True, bond_threshold=1.9)

or, with the pieces exposed::

    from pdb2pov import ParseOptions, SceneOptions, read_structure
    from pdb2pov import find_bonds, prepare_structure, write_scene

    structure, stats = read_structure("1crn.cif.gz", ParseOptions(chains="A"))
    options = SceneOptions(ball_stick=True, name="crambin")
    prepare_structure(structure, options)
    bonds = find_bonds(structure, options.bond_threshold)
    write_scene(structure, options, "crambin.pov", bonds)

Copyright (c) 1993-2026 Eric G. Suchanek, Ph.D.  Subject to the GNU License.
"""

from __future__ import annotations

import os

from .bonds import find_bonds
from .elements import (
    ELEMENT_UNKNOWN,
    ELEMENTS,
    element_count,
    element_index,
    element_pov_suffix,
    element_symbol,
    infer_element,
)
from .geometry import rotate_structure
from .mmcif import read_mmcif
from .options import (
    ALL_MODELS,
    AltLocPolicy,
    Backdrop,
    BondMode,
    Ground,
    InputFormat,
    ParseOptions,
    RadiiSet,
    SceneOptions,
)
from .readers import (
    ParseError,
    assign_types,
    read_atm,
    read_pdb,
    read_structure,
)
from .scene import (
    PDB2POV_VERSION,
    POV_VERSION,
    pov_identifier,
    prepare_structure,
    scene_text,
    write_scene,
)
from .structure import Atom, Extents, ParseStats, Structure

__version__ = "2.2.0"

__all__ = [
    "__version__",
    "PDB2POV_VERSION",
    "POV_VERSION",
    "ALL_MODELS",
    "ELEMENTS",
    "ELEMENT_UNKNOWN",
    "AltLocPolicy",
    "Atom",
    "Backdrop",
    "BondMode",
    "Extents",
    "Ground",
    "InputFormat",
    "ParseError",
    "ParseOptions",
    "ParseStats",
    "RadiiSet",
    "SceneOptions",
    "Structure",
    "assign_types",
    "convert",
    "element_count",
    "element_index",
    "element_pov_suffix",
    "element_symbol",
    "find_bonds",
    "include_dir",
    "infer_element",
    "pov_identifier",
    "prepare_structure",
    "read_atm",
    "read_mmcif",
    "read_pdb",
    "read_structure",
    "rotate_structure",
    "scene_text",
    "write_scene",
]


def include_dir() -> str:
    """
    Directory holding the bundled POV-Ray include files.

    The scenes reference ``atoms2.inc`` and one of the radius sets by name,
    so POV-Ray needs this on its library path::

        povray +Icrambin.pov +L$(python -m pdb2pov --include-dir)

    Shipping them inside the package is what makes this directory
    self-contained: copy it anywhere and the scenes it writes still render.
    """
    return os.path.join(os.path.dirname(os.path.abspath(__file__)), "include")


def convert(
    input_path: str,
    output_stem: str,
    parse_options: ParseOptions | None = None,
    scene_options: SceneOptions | None = None,
    **overrides,
) -> tuple[Structure, ParseStats, str]:
    """
    Read a structure and write a scene, the whole pipeline in one call.

    Keyword overrides are applied to whichever options object declares them,
    so ``convert("1crn", "out", ball_stick=True, chains="A")`` works without
    building either object by hand.  Returns the prepared structure, the
    parse statistics, and the path written.
    """
    popt = parse_options or ParseOptions()
    sopt = scene_options or SceneOptions()

    for key, value in overrides.items():
        if hasattr(popt, key):
            setattr(popt, key, value)
        elif hasattr(sopt, key):
            setattr(sopt, key, value)
        else:
            raise TypeError(f"convert() got an unexpected keyword argument {key!r}")

    if scene_options is None and "name" not in overrides:
        sopt.name = pov_identifier(output_stem)
    if popt.legacy_elements:
        sopt.legacy_drop_unknown = True

    structure, stats = read_structure(input_path, popt)
    prepare_structure(structure, sopt)

    bonds = None
    if sopt.ball_stick:
        bonds = find_bonds(
            structure,
            sopt.bond_threshold,
            sopt.bond_mode,
            sopt.bond_tolerance,
        )

    path = output_stem
    if not path.lower().endswith((".pov", ".inc")):
        path += sopt.extension

    write_scene(structure, sopt, path, bonds)
    return structure, stats, path
