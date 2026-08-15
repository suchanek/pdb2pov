"""
The parsed structure: atoms, and what a parse had to skip or guess at.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field

from .elements import ELEMENT_UNKNOWN, element_symbol

__all__ = ["Atom", "Structure", "Extents", "ParseStats", "MAX_RAD", "WATER_RESIDUES"]

#: Atomic radius in angstroms used to pad the bounding box and the enclosing
#: sphere.  Carried over from the C, where the comment is longer: it is no
#: longer an upper bound on every radius (potassium is 2.75 A van der Waals),
#: but raising it would move the camera and the reported enclosing radius of
#: every existing scene for a quarter of an angstrom of framing margin.
MAX_RAD = 2.5

#: Residue names that mean water.  One definition, so the ``--no-water``
#: filter and :attr:`Atom.is_water` cannot disagree about what a water is.
WATER_RESIDUES = frozenset(("HOH", "WAT", "DOD", "H2O", "TIP3", "SOL"))


@dataclass(slots=True)
class Atom:
    """
    One coordinate record.

    Everything the readers can recover is kept, not only what the scene
    writer needs -- occupancy and B-factor cost nothing to carry and callers
    using this as a PDB reader want them.
    """

    x: float
    y: float
    z: float
    name: str = ""  #: atom name, trimmed
    element: str = ""  #: upper-case symbol; "" if the file did not say
    res_name: str = ""
    chain_id: str = " "
    res_seq: int = 0
    i_code: str = " "
    alt_loc: str = " "
    serial: int = 0
    occupancy: float = 1.0
    temp_factor: float = 0.0
    charge: float = 0.0
    hetatm: bool = False
    model: int = 1

    #: Index into :data:`pdb2pov.elements.ELEMENTS`, assigned by
    #: :func:`pdb2pov.readers.assign_types`.  :data:`ELEMENT_UNKNOWN` means
    #: the atom draws as the generic ``Atom_X``.
    type_index: int = ELEMENT_UNKNOWN

    #: True when :attr:`element` was inferred from the atom name rather than
    #: read from the file.  Reported, so a questionable scene can be traced
    #: back to a file with no element column.
    element_inferred: bool = False

    @property
    def position(self) -> tuple[float, float, float]:
        return (self.x, self.y, self.z)

    @property
    def is_hydrogen(self) -> bool:
        """
        Hydrogen by assigned element, not by name.

        The bonding rules cap hydrogens at one partner.  Testing the name
        instead -- which is what the C did before 2.1 -- capped mercury,
        hafnium and helium too, so a mercury bridging two cysteines came out
        with one bond.
        """
        return self.type_index != ELEMENT_UNKNOWN and element_symbol(self.type_index) == "H"

    @property
    def is_water(self) -> bool:
        return self.res_name.strip().upper() in WATER_RESIDUES


@dataclass
class Extents:
    """An axis-aligned bounding box, already padded by :data:`MAX_RAD`."""

    xmin: float = 0.0
    xmax: float = 0.0
    ymin: float = 0.0
    ymax: float = 0.0
    zmin: float = 0.0
    zmax: float = 0.0


@dataclass
class ParseStats:
    """
    What a parse discarded or had to guess at, so it can be reported rather
    than swallowed.  The 1993 code dropped atoms silently in three places.
    """

    accepted: int = 0
    skipped_altloc: int = 0
    skipped_chain: int = 0
    skipped_model: int = 0
    skipped_hetatm: int = 0
    skipped_water: int = 0
    skipped_malformed: int = 0
    freeform_fallback: int = 0
    no_element_column: int = 0
    inferred_elements: int = 0
    generic: int = 0

    #: Element symbols (or atom names, when there was no symbol) that fell
    #: through to the generic texture, de-duplicated and capped for display.
    generic_symbols: list[str] = field(default_factory=list)

    #: ``(line number, text)`` for the first few unparseable records.  The C
    #: reported only a count; a line number is the difference between "some
    #: records are broken" and knowing which.
    malformed_examples: list[tuple[int, str]] = field(default_factory=list)

    #: Models present in the file, in the order encountered.
    models_seen: list[int] = field(default_factory=list)

    #: ``(atom name, chosen, rejected)`` for names that could have been read
    #: as a two-letter element and were not.  See
    #: :func:`pdb2pov.elements.infer_element_verbose`.
    ambiguous_names: list[tuple[str, str, str]] = field(default_factory=list)

    max_reported_symbols: int = 12
    max_reported_examples: int = 5

    def note_generic(self, symbol: str) -> None:
        symbol = symbol.strip() or "?"
        if symbol not in self.generic_symbols:
            if len(self.generic_symbols) < self.max_reported_symbols:
                self.generic_symbols.append(symbol)

    def note_ambiguous(self, name: str, chosen: str, rejected: str) -> None:
        entry = (name, chosen, rejected)
        if entry not in self.ambiguous_names:
            if len(self.ambiguous_names) < self.max_reported_symbols:
                self.ambiguous_names.append(entry)

    def note_malformed(self, lineno: int, text: str) -> None:
        self.skipped_malformed += 1
        if len(self.malformed_examples) < self.max_reported_examples:
            self.malformed_examples.append((lineno, text.rstrip("\r\n")))

    @property
    def truncated_symbols(self) -> bool:
        return len(self.generic_symbols) >= self.max_reported_symbols

    def lines(self, legacy: bool = False) -> list[str]:
        """Human-readable summary lines; empty when there is nothing to say."""
        out: list[str] = []
        if self.skipped_altloc:
            out.append(
                f"  {self.skipped_altloc} atom(s) in alternate conformations skipped "
                "(--keep-altlocs to keep them)"
            )
        if self.skipped_chain:
            out.append(f"  {self.skipped_chain} atom(s) outside the requested chain(s) skipped")
        if self.skipped_model:
            out.append(f"  {self.skipped_model} atom(s) in other models skipped (--model to choose)")
        if self.skipped_hetatm:
            out.append(f"  {self.skipped_hetatm} HETATM record(s) skipped")
        if self.skipped_water:
            out.append(f"  {self.skipped_water} water atom(s) skipped")
        if self.skipped_malformed:
            out.append(f"  {self.skipped_malformed} unparseable coordinate record(s) skipped")
            for lineno, text in self.malformed_examples:
                out.append(f"      line {lineno}: {text[:72]}")
        if self.freeform_fallback:
            out.append(
                f"  {self.freeform_fallback} record(s) did not honour the PDB column "
                "layout; read free-form"
            )
        if self.no_element_column:
            out.append(
                f"  {self.no_element_column} record(s) have no element column; "
                f"{self.inferred_elements} inferred from the atom name"
            )
        if self.ambiguous_names:
            readings = ", ".join(
                f"'{name}' as {chosen} not {rejected}"
                for name, chosen, rejected in self.ambiguous_names
            )
            out.append(
                f"  {len(self.ambiguous_names)} atom name(s) read by the PDB column "
                f"rule where a two-letter element was also possible ({readings});"
            )
            out.append(
                "      add an element column in 77-78, or use --legacy-elements, "
                "to override"
            )
        if self.generic:
            what = "are dropped (--legacy-elements)" if legacy else "render as Atom_X"
            symbols = ", ".join(self.generic_symbols)
            if self.truncated_symbols:
                symbols += ", ..."
            out.append(f"  {self.generic} atom(s) have no dedicated texture and {what} ({symbols})")
        return out


@dataclass
class Structure:
    """A parsed structure and where it came from."""

    atoms: list[Atom] = field(default_factory=list)
    title: str = ""
    source: str = ""
    source_format: str = ""

    def __len__(self) -> int:
        return len(self.atoms)

    def __iter__(self):
        return iter(self.atoms)

    def __getitem__(self, index):
        return self.atoms[index]

    @property
    def natoms(self) -> int:
        return len(self.atoms)

    # -- geometry ------------------------------------------------------

    def extents(self) -> Extents:
        """
        Axis-aligned bounds padded by :data:`MAX_RAD` to allow for the atom
        spheres.  An empty structure gives a zero-size box at the origin
        rather than the C's 9999.9 sentinels leaking into the output.
        """
        if not self.atoms:
            return Extents()

        xs = [a.x for a in self.atoms]
        ys = [a.y for a in self.atoms]
        zs = [a.z for a in self.atoms]
        return Extents(
            xmin=min(xs) - MAX_RAD,
            xmax=max(xs) + MAX_RAD,
            ymin=min(ys) - MAX_RAD,
            ymax=max(ys) + MAX_RAD,
            zmin=min(zs) - MAX_RAD,
            zmax=max(zs) + MAX_RAD,
        )

    def enclosing_radius(self) -> float:
        """
        Radius of an origin-centred sphere containing every atom.  Assumes
        :meth:`center` has already run, as the C does.
        """
        if not self.atoms:
            return 0.0
        return max(math.sqrt(a.x * a.x + a.y * a.y + a.z * a.z) for a in self.atoms) + MAX_RAD

    def centroid(self) -> tuple[float, float, float]:
        if not self.atoms:
            return (0.0, 0.0, 0.0)
        n = float(len(self.atoms))
        return (
            sum(a.x for a in self.atoms) / n,
            sum(a.y for a in self.atoms) / n,
            sum(a.z for a in self.atoms) / n,
        )

    def center(self) -> tuple[float, float, float]:
        """Translate to the centroid, returning the shift applied."""
        cx, cy, cz = self.centroid()
        for a in self.atoms:
            a.x -= cx
            a.y -= cy
            a.z -= cz
        return (cx, cy, cz)

    def flip_z(self) -> None:
        """PDB coordinates are right-handed; POV-Ray's world is left-handed."""
        for a in self.atoms:
            a.z = -a.z

    # -- selection -----------------------------------------------------

    def filtered(self, predicate) -> "Structure":
        """A new structure holding the atoms for which ``predicate`` is true."""
        return Structure(
            atoms=[a for a in self.atoms if predicate(a)],
            title=self.title,
            source=self.source,
            source_format=self.source_format,
        )

    def chains(self) -> list[str]:
        seen: list[str] = []
        for a in self.atoms:
            if a.chain_id not in seen:
                seen.append(a.chain_id)
        return seen

    def element_counts(self) -> dict[str, int]:
        counts: dict[str, int] = {}
        for a in self.atoms:
            key = a.element or "?"
            counts[key] = counts.get(key, 0) + 1
        return dict(sorted(counts.items(), key=lambda kv: (-kv[1], kv[0])))
