"""
Element table, radii and element inference.

The table mirrors the ``ELEMENTS`` array in ``pdb2pov.c``: a PDB element
symbol paired with the POV-Ray identifier suffix used to build
``Atom_<suffix>`` and ``Atom_Glass_<suffix>``.  Order matters only in that
the ``make check`` element grid and the C program's type codes follow it, so
it is kept identical.

The radii are duplicated from the bundled include files.  POV-Ray reads them
from ``atoms_vdw.inc`` and friends at render time, so nothing in the emitted
scene depends on the copies here -- they exist so the covalent bonding mode
can size bonds without parsing a POV-Ray include, and so callers using the
library for something other than scene generation have them.  ``make check``
in ``python/Makefile`` re-parses the include files and asserts they agree.
"""

from __future__ import annotations

from typing import NamedTuple

__all__ = [
    "Element",
    "ELEMENTS",
    "ELEMENT_UNKNOWN",
    "element_count",
    "element_symbol",
    "element_pov_suffix",
    "element_index",
    "covalent_radius",
    "vdw_radius",
    "infer_element",
    "infer_element_verbose",
    "legacy_guess_symbol",
    "normalise_symbol",
]


class Element(NamedTuple):
    """One row of the drawable-element table."""

    symbol: str  #: PDB element symbol, upper case, as in columns 77-78
    pov_suffix: str  #: suffix of the POV-Ray identifier, e.g. "Zn"


#: Index returned for an element with no dedicated texture.  Such atoms are
#: drawn with the generic ``Atom_X`` sphere rather than being dropped.
ELEMENT_UNKNOWN = -1

ELEMENTS: tuple[Element, ...] = (
    # organic and biological
    Element("H", "H"),
    Element("C", "C"),
    Element("N", "N"),
    Element("O", "O"),
    Element("S", "S"),
    Element("P", "P"),
    Element("SE", "Se"),
    # halogens
    Element("F", "F"),
    Element("CL", "Cl"),
    Element("BR", "Br"),
    Element("I", "I"),
    # alkali and alkaline earth
    Element("LI", "Li"),
    Element("NA", "Na"),
    Element("K", "K"),
    Element("MG", "Mg"),
    Element("CA", "Ca"),
    # transition and heavy metals
    Element("MN", "Mn"),
    Element("FE", "Fe"),
    Element("CO", "Co"),
    Element("NI", "Ni"),
    Element("CU", "Cu"),
    Element("ZN", "Zn"),
    Element("MO", "Mo"),
    Element("W", "W"),
    Element("AG", "Ag"),
    Element("CD", "Cd"),
    Element("PT", "Pt"),
    Element("AU", "Au"),
    Element("HG", "Hg"),
    # other
    Element("B", "B"),
    Element("SI", "Si"),
    Element("AS", "As"),
    Element("XE", "Xe"),
)

_INDEX_BY_SYMBOL = {e.symbol: i for i, e in enumerate(ELEMENTS)}

#: Van der Waals radii, in angstroms -- the values in ``atoms_vdw.inc``.  The
#: original eight keep their 1994 figures (Pauling and the CAChe software);
#: the rest are Bondi's.
VDW_RADII: dict[str, float] = {
    "H": 1.2, "C": 1.7, "N": 1.54, "O": 1.4, "S": 1.8, "P": 1.9,
    "CA": 1.274, "FE": 2.0, "SE": 1.9, "F": 1.47, "CL": 1.75, "BR": 1.85,
    "I": 1.98, "LI": 1.82, "NA": 2.27, "K": 2.75, "MG": 1.73, "MN": 2.05,
    "CO": 2.0, "NI": 1.63, "CU": 1.4, "ZN": 1.39, "MO": 2.1, "W": 2.1,
    "AG": 1.72, "CD": 1.58, "PT": 1.75, "AU": 1.66, "HG": 1.55, "B": 1.92,
    "SI": 2.1, "AS": 1.85, "XE": 2.16,
}
VDW_GENERIC = 1.7

#: Covalent radii, in angstroms -- the values in ``atoms_covalent.inc``.
#: Cordero et al. (2008) for everything added in 2.2, low-spin where a metal
#: has spin-state-dependent radii.
COVALENT_RADII: dict[str, float] = {
    "H": 0.35, "C": 0.77, "N": 0.74, "O": 0.73, "S": 1.04, "P": 1.1,
    "CA": 1.74, "FE": 1.17, "SE": 1.2, "F": 0.57, "CL": 1.02, "BR": 1.2,
    "I": 1.39, "LI": 1.28, "NA": 1.66, "K": 2.03, "MG": 1.41, "MN": 1.61,
    "CO": 1.26, "NI": 1.24, "CU": 1.32, "ZN": 1.22, "MO": 1.54, "W": 1.62,
    "AG": 1.45, "CD": 1.44, "PT": 1.36, "AU": 1.36, "HG": 1.32, "B": 0.84,
    "SI": 1.11, "AS": 1.19, "XE": 1.4,
}
COVALENT_GENERIC = 1.0


def element_count() -> int:
    """Number of elements with a dedicated texture."""
    return len(ELEMENTS)


def element_symbol(index: int) -> str:
    """PDB symbol for a table index; ``"X"`` for :data:`ELEMENT_UNKNOWN`."""
    if 0 <= index < len(ELEMENTS):
        return ELEMENTS[index].symbol
    return "X"


def element_pov_suffix(index: int) -> str:
    """POV-Ray identifier suffix for a table index; ``"X"`` if unlisted."""
    if 0 <= index < len(ELEMENTS):
        return ELEMENTS[index].pov_suffix
    return "X"


def element_index(symbol: str | None) -> int:
    """
    Table index for a PDB element symbol, or :data:`ELEMENT_UNKNOWN`.

    Deuterium and tritium are hydrogen as far as a picture is concerned.
    """
    if not symbol:
        return ELEMENT_UNKNOWN

    sym = normalise_symbol(symbol)
    if sym in ("D", "T"):
        sym = "H"
    return _INDEX_BY_SYMBOL.get(sym, ELEMENT_UNKNOWN)


def normalise_symbol(symbol: str) -> str:
    """
    Reduce a raw element field to a bare upper-case symbol.

    mmCIF and some hand-made PDB files carry a charge in the same field
    (``FE2+``, ``CL-``), and PDB columns 77-78 are blank-padded.  Only the
    leading letters are the symbol.
    """
    out = []
    for ch in symbol.strip():
        if ch.isalpha():
            out.append(ch)
        else:
            break
    return "".join(out).upper()


def vdw_radius(symbol: str) -> float:
    """Van der Waals radius in angstroms, generic if the element is unlisted."""
    return VDW_RADII.get(normalise_symbol(symbol), VDW_GENERIC)


def covalent_radius(symbol: str) -> float:
    """Covalent radius in angstroms, generic if the element is unlisted."""
    return COVALENT_RADII.get(normalise_symbol(symbol), COVALENT_GENERIC)


# ----------------------------------------------------------------------
# Element inference
# ----------------------------------------------------------------------

#: Every element symbol, not just the drawable ones.  Inference has to be
#: able to recognise, say, gadolinium in order to decide that ``GD`` in
#: columns 13-14 is one atom of it rather than a carbon named "GD" -- even
#: though the result is a grey ``Atom_X``.
_ALL_SYMBOLS = frozenset(
    """
    H HE LI BE B C N O F NE NA MG AL SI P S CL AR K CA SC TI V CR MN FE CO NI
    CU ZN GA GE AS SE BR KR RB SR Y ZR NB MO TC RU RH PD AG CD IN SN SB TE I
    XE CS BA LA CE PR ND PM SM EU GD TB DY HO ER TM YB LU HF TA W RE OS IR PT
    AU HG TL PB BI PO AT RN FR RA AC TH PA U NP PU AM CM BK CF ES FM MD NO LR
    RF DB SG BH HS MT DS RG CN NH FL MC LV TS OG D T
    """.split()
)

_TWO_LETTER_SYMBOLS = frozenset(s for s in _ALL_SYMBOLS if len(s) == 2)

#: Residues whose atom names follow the protein/nucleic naming scheme, where
#: the element is always the first letter of the trimmed name.  Getting this
#: right is what stops a haem nitrogen named ``NA`` becoming sodium and a
#: glutamate delta carbon named ``CD`` becoming cadmium.
_STANDARD_RESIDUES = frozenset(
    """
    ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER THR TRP
    TYR VAL SEC PYL ASX GLX UNK
    A C G U I DA DC DG DT DI DU N
    HOH WAT DOD H2O
    """.split()
)

#: Selenomethionine and selenocysteine are standard enough to be named like
#: amino acids while still carrying a two-letter element.
_STANDARD_TWO_LETTER = {
    "MSE": ("SE",),
    "SEC": ("SE",),
    "CSE": ("SE",),
}


def infer_element(name_raw: str, res_name: str = "", hetatm: bool = False) -> str:
    """Element symbol inferred from an atom name; see :func:`infer_element_verbose`."""
    return infer_element_verbose(name_raw, res_name, hetatm)[0]


def infer_element_verbose(
    name_raw: str,
    res_name: str = "",
    hetatm: bool = False,
) -> tuple[str, str | None]:
    """
    Work out an element symbol from an atom name when no element column exists.

    ``name_raw`` is the *untrimmed* four-character atom-name field, because
    its alignment carries the answer: the PDB format right-justifies a
    two-letter element in columns 13-14, so a name starting in column 13 is
    the strong signal that the first two characters are the symbol.  Callers
    that only have a trimmed name should pass it as-is; inference then falls
    back to the weaker rules below.

    This replaces the pre-2.1 guess, which looked at the first one or two
    characters of the name and knew seven letters: sodium read as nitrogen,
    chlorine as carbon, fluorine as iron, and everything else was dropped.
    See :func:`legacy_guess_symbol`, kept for ``--legacy-elements``.

    Returns ``(symbol, rejected)``.  ``symbol`` is ``""`` when nothing can be
    concluded, which the caller reports rather than guessing at.  ``rejected``
    names the other reading when the name could plausibly have been a
    two-letter element and was not read as one -- ``(" FE ", "XXX")`` gives
    ``("F", "FE")`` -- so an unusual file gets a warning instead of a silent
    fluorine.
    """
    padded = f"{name_raw:<4}"[:4] if len(name_raw) < 4 else name_raw
    name = name_raw.strip().upper()
    res = res_name.strip().upper()

    if not name:
        return ("", None)

    # Hydrogen names are often digit-prefixed: 1HB, 2HG1, 3HD2.
    if name[0].isdigit():
        rest = name.lstrip("0123456789")
        return (rest[:1] if rest[:1] in _ALL_SYMBOLS else "", None)

    # A monatomic ion is written with the element as both residue and atom
    # name -- ZN/ZN, CU/CU, NA/NA.  That agreement is worth more than column
    # alignment, which hand-edited files get wrong in exactly this case.
    alpha = "".join(ch for ch in name if ch.isalpha())
    if alpha and alpha == res and alpha in _ALL_SYMBOLS:
        return (alpha, None)

    if res in _STANDARD_TWO_LETTER:
        for sym in _STANDARD_TWO_LETTER[res]:
            if alpha.startswith(sym):
                return (sym, None)

    # In a standard residue the element is the leading letter, full stop.
    if res in _STANDARD_RESIDUES:
        return (name[0] if name[0] in _ALL_SYMBOLS else "", None)

    # Otherwise trust the column layout: a two-letter element is
    # right-justified into columns 13-14, so it starts in column 13.
    if padded[0] not in " \t":
        two = padded[:2].strip().upper()
        if len(two) == 2 and two in _TWO_LETTER_SYMBOLS:
            return (two, None)

    # Name starts in column 14, which the format reserves for one-letter
    # elements.  Deliberately not second-guessed: " NA " inside a haem is
    # nitrogen, and preferring the two-letter reading here would break far
    # more structures than the hand-edited files it would fix.  It is
    # reported, though, because a hand-edited " FE " does exist.
    if alpha[:1] in _ALL_SYMBOLS:
        rejected = alpha[:2] if alpha[:2] in _TWO_LETTER_SYMBOLS else None
        return (alpha[0], rejected)

    return ("", None)


def legacy_guess_symbol(name: str, hetatm: bool = False) -> str:
    """
    The pre-2.1 element guess, reproduced exactly for ``--legacy-elements``.

    It knows seven first letters and is wrong for any two-letter element
    sharing a first letter with a one-letter one.  ``hetatm`` drives the
    alpha-carbon hack: in an ATOM record a name beginning ``CA`` is
    C-alpha, not calcium, which the C code arranged by blanking the second
    character of the name before guessing.
    """
    if not name:
        return ""

    upper = name.upper()
    a = upper[0]
    b = upper[1] if len(upper) > 1 else " "

    if not hetatm and upper.startswith("CA"):
        b = " "  # the alpha-carbon hack

    if a == "H":
        return "H"
    if a == "C":
        return "CA" if b == "A" else "C"
    if a == "O":
        return "O"
    if a == "N":
        return "N"
    if a == "S":
        return "S"
    if a == "P":
        return "P"
    if a == "F":
        return "FE"
    return ""
