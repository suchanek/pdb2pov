"""
Bond search.

Two things differ from the C, neither of them visible in the output.

*Cost.*  The C compares every pair, twice -- once to count and once to
record.  That is fine for crambin's 327 atoms and hopeless for a ribosome:
100,000 atoms is five billion pair tests.  Here the atoms go into a uniform
grid of cells one cutoff wide, so each atom only looks at the 27 cells around
it and the search is linear in the number of atoms.

*Order.*  The result is nevertheless identical, pair for pair and in the same
sequence, because the traversal order the C's hydrogen rule depends on is
reproduced deliberately: the outer index runs down from the last atom, the
inner index runs down from just below it, and a hydrogen stops at its first
partner.  Hydrogens conventionally follow the heavy atom they hang off, so
working backwards finds that heavy atom first.

:func:`find_bonds` also gains a covalent mode, in which the cutoff is the sum
of the two covalent radii rather than one distance for the whole structure.
"""

from __future__ import annotations

import math
from collections import defaultdict

from .elements import ELEMENT_UNKNOWN, covalent_radius, element_symbol
from .options import BondMode
from .structure import Structure

__all__ = ["find_bonds", "Bond"]

Bond = tuple[int, int]


def _build_grid(structure: Structure, cell: float) -> dict[tuple[int, int, int], list[int]]:
    grid: dict[tuple[int, int, int], list[int]] = defaultdict(list)
    inv = 1.0 / cell
    for index, atom in enumerate(structure.atoms):
        key = (
            math.floor(atom.x * inv),
            math.floor(atom.y * inv),
            math.floor(atom.z * inv),
        )
        grid[key].append(index)
    return grid


_NEIGHBOUR_OFFSETS = tuple(
    (dx, dy, dz)
    for dx in (-1, 0, 1)
    for dy in (-1, 0, 1)
    for dz in (-1, 0, 1)
)


def _bonding_symbol(atom) -> str:
    """The symbol a bond radius should be looked up under."""
    if atom.element:
        return atom.element
    if atom.type_index != ELEMENT_UNKNOWN:
        return element_symbol(atom.type_index)
    return ""


def find_bonds(
    structure: Structure,
    threshold: float | None = 2.2,
    mode: BondMode = BondMode.DISTANCE,
    tolerance: float = 0.45,
    max_bonds: int | None = None,
) -> list[Bond]:
    """
    Bonded pairs as ``(lower index, higher index)``.

    ``threshold`` is the cutoff in angstroms for :data:`BondMode.DISTANCE`,
    where it is required.  In :data:`BondMode.COVALENT` a pair bonds when its
    separation is at most the sum of the two covalent radii plus
    ``tolerance``, and ``threshold`` becomes an optional upper bound on the
    search: pass ``None`` to let the radii decide entirely.
    """
    atoms = structure.atoms
    n = len(atoms)
    if n < 2:
        return []

    covalent = mode is BondMode.COVALENT
    if covalent:
        radii = [covalent_radius(_bonding_symbol(a)) for a in atoms]
        search = 2.0 * max(radii) + tolerance
        if threshold is not None:
            search = min(search, threshold)
    else:
        if threshold is None:
            raise ValueError("a bond threshold is required in distance mode")
        radii = []
        search = threshold

    if search <= 0.0:
        return []

    grid = _build_grid(structure, search)
    inv = 1.0 / search
    is_h = [a.is_hydrogen for a in atoms]

    bonds: list[Bond] = []

    for i in range(n - 1, 0, -1):
        ai = atoms[i]
        cx = math.floor(ai.x * inv)
        cy = math.floor(ai.y * inv)
        cz = math.floor(ai.z * inv)

        candidates: list[int] = []
        for dx, dy, dz in _NEIGHBOUR_OFFSETS:
            cellmates = grid.get((cx + dx, cy + dy, cz + dz))
            if cellmates:
                candidates.extend(j for j in cellmates if j < i)

        if not candidates:
            continue

        # Descending, so the traversal matches the C's inner loop and the
        # hydrogen rule below picks the same partner it always has.
        candidates.sort(reverse=True)

        hi = is_h[i]
        for j in candidates:
            if hi and is_h[j]:
                continue  # never bond hydrogens to each other

            aj = atoms[j]

            # The same sequence of rejections the C uses, hypot included:
            # at the cutoff a differently associated sqrt can land on the
            # other side of the comparison, and these scenes are diffed
            # against the C's byte for byte.
            dx = abs(ai.x - aj.x)
            if dx > search:
                continue
            dy = abs(ai.y - aj.y)
            if dy > search:
                continue
            dxy = math.hypot(dx, dy)
            if dxy > search:
                continue
            dz = abs(ai.z - aj.z)
            if dz > search:
                continue
            dist = math.hypot(dxy, dz)
            if dist > search:
                continue

            if covalent and dist > radii[i] + radii[j] + tolerance:
                continue

            bonds.append((j, i))
            if max_bonds is not None and len(bonds) >= max_bonds:
                return bonds

            if hi:
                break  # one bond per hydrogen

    return bonds
