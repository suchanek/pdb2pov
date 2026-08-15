"""The bond search: same answers as the C, in the same order, far faster."""

from __future__ import annotations

import math
import random

from pdb2pov import BondMode, ParseOptions, find_bonds, read_structure
from pdb2pov.structure import Atom, Structure


def brute_force(structure: Structure, threshold: float) -> list[tuple[int, int]]:
    """The C's loop, transcribed, as the reference the grid must agree with."""
    bonds = []
    atoms = structure.atoms
    for i in range(len(atoms) - 1, 0, -1):
        for j in range(i - 1, -1, -1):
            if atoms[i].is_hydrogen and atoms[j].is_hydrogen:
                continue
            dx = abs(atoms[i].x - atoms[j].x)
            if dx > threshold:
                continue
            dy = abs(atoms[i].y - atoms[j].y)
            if dy > threshold:
                continue
            dxy = math.hypot(dx, dy)
            if dxy > threshold:
                continue
            dz = abs(atoms[i].z - atoms[j].z)
            if dz > threshold:
                continue
            if math.hypot(dxy, dz) > threshold:
                continue
            bonds.append((j, i))
            if atoms[i].is_hydrogen:
                break
    return bonds


def random_structure(n: int, seed: int, box: float = 12.0) -> Structure:
    rng = random.Random(seed)
    atoms = []
    for index in range(n):
        element = "H" if index % 4 == 0 else "C"
        atoms.append(
            Atom(
                x=rng.uniform(-box, box),
                y=rng.uniform(-box, box),
                z=rng.uniform(-box, box),
                name=element,
                element=element,
            )
        )
    structure = Structure(atoms=atoms)
    from pdb2pov import assign_types

    assign_types(structure)
    return structure


def test_grid_search_matches_the_brute_force_loop_exactly():
    for seed in range(6):
        structure = random_structure(220, seed)
        for threshold in (1.2, 2.2, 3.7):
            assert find_bonds(structure, threshold) == brute_force(structure, threshold)


def test_the_first_atom_can_be_bonded():
    """The C's inner loop ran `j > 0` until 2.0, so atom 0 never bonded."""
    structure = Structure(
        atoms=[
            Atom(0.0, 0.0, 0.0, name="N", element="N"),
            Atom(1.4, 0.0, 0.0, name="CA", element="C"),
        ]
    )
    from pdb2pov import assign_types

    assign_types(structure)
    assert find_bonds(structure, 2.2) == [(0, 1)]


def test_hydrogens_take_one_partner_and_never_each_other():
    from pdb2pov import assign_types

    structure = Structure(
        atoms=[
            Atom(0.0, 0.0, 0.0, name="C", element="C"),
            Atom(1.1, 0.0, 0.0, name="H1", element="H"),
            Atom(1.1, 0.4, 0.0, name="H2", element="H"),
        ]
    )
    assign_types(structure)

    bonds = find_bonds(structure, 2.2)
    assert bonds == [(0, 2), (0, 1)]  # no H-H, one bond each


def test_mercury_is_not_capped_like_a_hydrogen():
    """Hydrogen is decided by element, not by a name beginning with H."""
    from pdb2pov import assign_types

    structure = Structure(
        atoms=[
            Atom(-2.0, 0.0, 0.0, name="SG", element="S"),
            Atom(2.0, 0.0, 0.0, name="SG", element="S"),
            Atom(0.0, 0.0, 0.0, name="HG", element="HG"),
        ]
    )
    assign_types(structure)

    assert len(find_bonds(structure, 2.2)) == 2


def test_covalent_mode_finds_long_and_short_bonds_at_once():
    """One cutoff cannot cover a 1.1 A C-H and a 2.05 A S-S; radii can."""
    from pdb2pov import assign_types

    structure = Structure(
        atoms=[
            Atom(0.0, 0.0, 0.0, name="C", element="C"),
            Atom(1.1, 0.0, 0.0, name="H", element="H"),
            Atom(0.0, 2.05, 0.0, name="SG", element="S"),
            Atom(0.0, 4.10, 0.0, name="SG", element="S"),
        ]
    )
    assign_types(structure)

    covalent = find_bonds(structure, None, BondMode.COVALENT)
    pairs = {tuple(sorted(b)) for b in covalent}

    assert (0, 1) in pairs  # C-H
    assert (2, 3) in pairs  # S-S
    assert (1, 2) not in pairs  # H to a sulphur 2.3 A away is not a bond


def test_crambin_bond_count_is_unchanged(crambin_pdb):
    from pdb2pov import prepare_structure
    from pdb2pov.options import SceneOptions

    structure, _ = read_structure(crambin_pdb, ParseOptions())
    prepare_structure(structure, SceneOptions())

    assert find_bonds(structure, 1.9) == brute_force(structure, 1.9)


def test_the_search_is_linear_enough_to_handle_a_large_structure():
    """
    50,000 atoms is a ribosome-sized job the C would take five hours over.
    This is a smoke test, not a benchmark: it only has to finish.
    """
    structure = random_structure(50_000, seed=99, box=90.0)
    bonds = find_bonds(structure, 2.2)
    assert bonds  # and, mostly, that we got here at all
