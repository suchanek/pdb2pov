"""
Readers for PDB, PDBx/mmCIF and the old ``.atm`` format.

What these do that ``pdb2pov.c`` does not:

* **mmCIF.**  Structures with more than 99,999 atoms or more than 62 chains
  are only distributed as mmCIF; the PDB format cannot express them.  See
  :mod:`pdb2pov.mmcif`.
* **Compression.**  ``.gz``, ``.bz2`` and ``.xz`` are opened transparently,
  which is how the wwPDB actually ships files.
* **Models.**  ``MODEL``/``ENDMDL`` are understood, so an NMR ensemble
  converts one model at a time instead of by accident.  The C stopped at the
  first line beginning ``END``, which happens to be ``ENDMDL`` -- right
  answer, wrong reason, and no way to ask for model 7.
* **No fixed limits.**  The C read into a 256-byte line buffer and a
  pre-counted fixed array; a longer line was silently split into two records.
* **Element inference** that knows about column alignment and residue names
  rather than seven first letters.  See :func:`pdb2pov.elements.infer_element`.
* **Diagnostics with line numbers**, and a strict mode that refuses to
  quietly drop records.

Everything the C accepted is still accepted, including its free-form
coordinate fallback and the non-standard trailing charge column.
"""

from __future__ import annotations

import bz2
import gzip
import lzma
import os
import sys
from contextlib import contextmanager
from typing import IO, Iterable, Iterator

from .elements import (
    ELEMENT_UNKNOWN,
    element_index,
    infer_element,
    infer_element_verbose,
    legacy_guess_symbol,
    normalise_symbol,
)
from .options import ALL_MODELS, AltLocPolicy, InputFormat, ParseOptions
from .structure import WATER_RESIDUES, Atom, ParseStats, Structure

__all__ = [
    "ParseError",
    "read_structure",
    "read_pdb",
    "read_atm",
    "assign_types",
    "apply_altloc_policy",
    "detect_format",
    "resolve_input_path",
    "open_text",
    "decode_hybrid36",
]


class ParseError(ValueError):
    """Raised in strict mode when a coordinate record will not parse."""


# ----------------------------------------------------------------------
# Opening files
# ----------------------------------------------------------------------

_COMPRESSORS = {".gz": gzip.open, ".bz2": bz2.open, ".xz": lzma.open}

#: Extensions tried, in order, when the path given has none that exists.
#: ``.pdb`` first keeps the C's behaviour of appending it to a bare stem.
_INPUT_SUFFIXES = (
    ".pdb", ".ent", ".cif", ".mmcif", ".atm",
    ".pdb.gz", ".ent.gz", ".cif.gz", ".mmcif.gz",
)


def resolve_input_path(path: str, fmt: InputFormat = InputFormat.AUTO) -> str:
    """
    Find the file the user meant.

    ``pdb2pov`` has always taken a stem and appended an extension.  That is
    kept -- ``pdb2pov 1crn out`` still reads ``1crn.pdb`` -- but a path that
    exists as given is used as given, so ``1crn.cif.gz`` works too.  ``-``
    means standard input.
    """
    if path == "-":
        return path
    if os.path.exists(path):
        return path

    preferred: tuple[str, ...] = ()
    if fmt is InputFormat.ATM:
        preferred = (".atm",)
    elif fmt is InputFormat.CIF:
        preferred = (".cif", ".mmcif", ".cif.gz", ".mmcif.gz")
    elif fmt is InputFormat.PDB:
        preferred = (".pdb", ".ent", ".pdb.gz", ".ent.gz")

    for suffix in preferred + _INPUT_SUFFIXES:
        candidate = path + suffix
        if os.path.exists(candidate):
            return candidate

    # Nothing matched; report the historical guess so the error names the
    # file the user probably meant.
    default = ".atm" if fmt is InputFormat.ATM else ".pdb"
    return path + default


@contextmanager
def open_text(path: str) -> Iterator[IO[str]]:
    """
    Open a possibly compressed structure file for reading text.

    Decoding is UTF-8 with replacement.  Coordinate files are ASCII, but
    depositor names in the header are not always, and a stray byte in a
    REMARK is no reason to refuse to draw the molecule.
    """
    if path == "-":
        yield sys.stdin
        return

    for suffix, opener in _COMPRESSORS.items():
        if path.endswith(suffix):
            with opener(path, "rt", encoding="utf-8", errors="replace") as handle:
                yield handle
            return

    with open(path, "r", encoding="utf-8", errors="replace") as handle:
        yield handle


def detect_format(path: str, sample: Iterable[str] = ()) -> InputFormat:
    """
    Decide what a file is, by extension first and content second.

    Content wins over a misleading name: a ``.pdb`` holding mmCIF is read as
    mmCIF rather than yielding "no atoms found".
    """
    lowered = path.lower()
    for suffix in _COMPRESSORS:
        if lowered.endswith(suffix):
            lowered = lowered[: -len(suffix)]

    guess = InputFormat.AUTO
    if lowered.endswith((".cif", ".mmcif")):
        guess = InputFormat.CIF
    elif lowered.endswith(".atm"):
        guess = InputFormat.ATM
    elif lowered.endswith((".pdb", ".ent")):
        guess = InputFormat.PDB

    for line in sample:
        stripped = line.strip()
        if not stripped:
            continue
        if stripped.startswith(("data_", "loop_", "_atom_site.", "#")):
            return InputFormat.CIF
        if line.startswith(("ATOM", "HETATM", "HEADER", "REMARK", "MODEL", "CRYST1", "TITLE")):
            return InputFormat.PDB

    return guess if guess is not InputFormat.AUTO else InputFormat.PDB


# ----------------------------------------------------------------------
# Small field helpers
# ----------------------------------------------------------------------

_HYBRID36_UPPER_OFFSET = {w: 10 * 36 ** (w - 1) - 10 ** w for w in range(2, 8)}


def decode_hybrid36(text: str, width: int) -> int:
    """
    Decode a PDB serial or residue number, hybrid-36 included.

    Files with more than 99,999 atoms encode the overflow in base 36 --
    ``A0000`` is 100,000.  Plain ``int()`` raises on those, which is enough
    to make a large structure look malformed.
    """
    token = text.strip()
    if not token:
        return 0
    if token[0].isdigit() or token[0] in "+-":
        try:
            return int(token)
        except ValueError:
            return 0

    offset = _HYBRID36_UPPER_OFFSET.get(width)
    if offset is None or not token.isalnum():
        return 0
    try:
        value = int(token, 36)
    except ValueError:
        return 0
    if token[0].isupper():
        return value - offset
    return value - offset + 26 * 36 ** (width - 1)


def _float_or(text: str, default: float) -> float:
    token = text.strip()
    if not token:
        return default
    try:
        return float(token)
    except ValueError:
        return default


def _formal_charge(text: str) -> float:
    """
    PDB columns 79-80 hold a formal charge as ``2+`` or ``1-``.  mmCIF spells
    the same thing ``2`` or ``-1``.  Both are accepted.
    """
    token = text.strip()
    if not token:
        return 0.0
    if token[-1] in "+-":
        token = token[-1] + token[:-1]
    try:
        return float(token)
    except ValueError:
        return 0.0


# ----------------------------------------------------------------------
# PDB
# ----------------------------------------------------------------------

def _is_coordinate_record(line: str) -> bool:
    upper = line[:6].upper()
    return upper.startswith("ATOM") or upper.startswith("HETATM")


def read_pdb(
    lines: Iterable[str],
    options: ParseOptions | None = None,
    stats: ParseStats | None = None,
) -> tuple[Structure, ParseStats]:
    """
    Read PDB coordinate records from an iterable of lines.

    Sample input, all of which the C accepted and which still parses here::

        ATOM     21  N   CYS     3      13.507  11.238   8.398
        HETATM   99  FE  XXX     4      33.265  10.443  10.307
        HETATM       CU+2CU    152     -20.577   2.601 -14.587       2.000

    The trailing charge column in the third record is not genuine PDB and is
    optional.
    """
    opt = options or ParseOptions()
    st = stats or ParseStats()
    structure = Structure(source_format="pdb")
    atoms = structure.atoms

    model = 1
    first_model: int | None = None
    want_all_models = opt.model == ALL_MODELS

    for lineno, raw in enumerate(lines, start=1):
        line = raw.rstrip("\n").rstrip("\r")
        record = line[:6].upper()

        if record.startswith("END") and not record.startswith("ENDMDL"):
            break

        if record.startswith("MODEL"):
            model = decode_hybrid36(line[10:14], 4) or (model + 1)
            if model not in st.models_seen:
                st.models_seen.append(model)
            if first_model is None:
                first_model = model
            continue

        if record.startswith("ENDMDL"):
            continue

        if line.strip() == ".":
            break

        if not _is_coordinate_record(line):
            continue

        if first_model is None:
            first_model = model
            if model not in st.models_seen:
                st.models_seen.append(model)

        wanted_model = opt.model if opt.model is not None else first_model
        if not want_all_models and model != wanted_model:
            st.skipped_model += 1
            continue

        hetatm = record.startswith("HETATM")

        # --- filters that need no parsing -----------------------------
        alt_loc = line[16] if len(line) > 16 else " "
        if opt.altloc is AltLocPolicy.A and alt_loc not in (" ", "A", "a", ""):
            st.skipped_altloc += 1
            continue

        chain_id = line[21] if len(line) > 21 else " "
        if not opt.accepts_chain(chain_id):
            st.skipped_chain += 1
            continue

        if hetatm and not opt.keep_hetatm:
            st.skipped_hetatm += 1
            continue

        res_name = line[17:20].strip()
        if not opt.keep_water and res_name.upper() in WATER_RESIDUES:
            st.skipped_water += 1
            continue

        # --- the record proper ----------------------------------------
        name_field = line[12:16]
        name = name_field.strip()
        if not name:
            _malformed(st, opt, lineno, line, "no atom name in columns 13-16")
            continue

        coords = _read_coordinates(line)
        if coords is None:
            _malformed(st, opt, lineno, line, "no coordinates in columns 31-54")
            continue
        (x, y, z), freeform = coords
        if freeform:
            st.freeform_fallback += 1

        element = normalise_symbol(line[76:78]) if len(line) > 76 else ""
        inferred = False
        if not element:
            st.no_element_column += 1
            element, rejected = infer_element_verbose(name_field, res_name, hetatm)
            if element:
                inferred = True
                st.inferred_elements += 1
            if rejected:
                st.note_ambiguous(name, element, rejected)

        charge = _formal_charge(line[78:80]) if len(line) > 78 else 0.0
        if charge == 0.0 and len(line) > 70:
            # The non-standard trailing charge the 1993 parser read.  Only
            # consulted when the standard column said nothing, and only when
            # it does not collide with the element column.
            legacy_charge = line[70:76].strip()
            if legacy_charge:
                charge = _float_or(legacy_charge, 0.0)

        atoms.append(
            Atom(
                x=x,
                y=y,
                z=z,
                name=name,
                element=element,
                res_name=res_name,
                chain_id=chain_id,
                res_seq=decode_hybrid36(line[22:26], 4),
                i_code=line[26] if len(line) > 26 else " ",
                alt_loc=alt_loc or " ",
                serial=decode_hybrid36(line[6:11], 5),
                occupancy=_float_or(line[54:60], 1.0),
                temp_factor=_float_or(line[60:66], 0.0),
                charge=charge,
                hetatm=hetatm,
                model=model,
                element_inferred=inferred,
            )
        )

        if opt.max_atoms is not None and len(atoms) >= opt.max_atoms:
            break

    apply_altloc_policy(structure, opt, st)
    st.accepted = len(structure.atoms)
    return structure, st


def _read_coordinates(line: str) -> tuple[tuple[float, float, float], bool] | None:
    """
    Coordinates from columns 31-54, or free-form from column 31 onwards.

    Fixed columns are what the format specifies.  Hand-edited files exist
    that do not honour them, and the C fell back to a free-form scan for
    those, so nothing that converted before stops converting now.  The
    fallback is reported.
    """
    try:
        return (
            (float(line[30:38]), float(line[38:46]), float(line[46:54])),
            False,
        )
    except (ValueError, IndexError):
        pass

    tail = line[30:].split()
    if len(tail) >= 3:
        try:
            return ((float(tail[0]), float(tail[1]), float(tail[2])), True)
        except ValueError:
            return None
    return None


def _malformed(
    st: ParseStats, opt: ParseOptions, lineno: int, line: str, why: str
) -> None:
    if opt.strict:
        raise ParseError(f"line {lineno}: {why}: {line[:72]!r}")
    st.note_malformed(lineno, line)


def apply_altloc_policy(
    structure: Structure, opt: ParseOptions, st: ParseStats
) -> None:
    """
    Resolve alternate conformations for the policies that need every
    conformer in hand: keep the first seen, or the highest occupancy.

    The choice is made **per residue**, not per atom, and the residue's name
    is deliberately not part of the grouping key.  Microheterogeneity is why:
    in 1CBN residue 22 is modelled as both proline and serine and residue 25
    as both leucine and isoleucine, sharing one sequence position.  Choosing
    atom by atom would take proline's ring and serine's hydroxyl from the
    same position and draw a residue that does not exist.  Picking one
    altLoc letter for the whole residue cannot do that.

    ``A`` and ``ALL`` are decided record by record while reading, so they do
    not reach here.
    """
    if opt.altloc in (AltLocPolicy.A, AltLocPolicy.ALL):
        return

    def residue_of(atom: Atom) -> tuple:
        return (atom.model, atom.chain_id, atom.res_seq, atom.i_code)

    # residue -> letter -> [order first seen, occupancy total, atom count]
    conformers: dict[tuple, dict[str, list[float]]] = {}
    for atom in structure.atoms:
        if atom.alt_loc in (" ", ""):
            continue
        letters = conformers.setdefault(residue_of(atom), {})
        entry = letters.get(atom.alt_loc)
        if entry is None:
            letters[atom.alt_loc] = [float(len(letters)), atom.occupancy, 1.0]
        else:
            entry[1] += atom.occupancy
            entry[2] += 1.0

    if not conformers:
        return

    winners: dict[tuple, str] = {}
    for residue, letters in conformers.items():
        if opt.altloc is AltLocPolicy.FIRST:
            winners[residue] = min(letters, key=lambda ltr: letters[ltr][0])
        else:
            # Mean occupancy, not the total: conformers of a microheterogeneous
            # residue have different atom counts, and the bigger side-chain
            # should not win on size alone.  Ties go to the first seen.
            winners[residue] = max(
                letters,
                key=lambda ltr: (letters[ltr][1] / letters[ltr][2], -letters[ltr][0]),
            )

    kept: list[Atom] = []
    dropped = 0
    for atom in structure.atoms:
        if atom.alt_loc in (" ", "") or winners.get(residue_of(atom)) == atom.alt_loc:
            kept.append(atom)
        else:
            dropped += 1

    if dropped:
        st.skipped_altloc += dropped
        structure.atoms[:] = kept


# ----------------------------------------------------------------------
# .atm
# ----------------------------------------------------------------------


def read_atm(
    lines: Iterable[str],
    options: ParseOptions | None = None,
    stats: ParseStats | None = None,
) -> tuple[Structure, ParseStats]:
    """
    Read the alternative ``.atm`` format::

        1 0x7042  Ca  0.0335223  5.4441102  0.0069856   20    0   sp   Ca  0x3

    The format has no element column, so types come from the atom name.
    """
    opt = options or ParseOptions()
    st = stats or ParseStats()
    structure = Structure(source_format="atm")

    for lineno, raw in enumerate(lines, start=1):
        line = raw.rstrip("\n").rstrip("\r")
        if line[:3].upper() == "END" or line.strip() == ".":
            break
        if not line.strip():
            continue

        fields = line.split()
        if len(fields) < 6:
            _malformed(st, opt, lineno, line, "expected at least six fields")
            continue

        try:
            serial = int(fields[0], 0)
        except ValueError:
            _malformed(st, opt, lineno, line, "first field is not an atom number")
            continue

        name = fields[2][:4]
        try:
            x, y, z = (float(fields[3]), float(fields[4]), float(fields[5]))
        except ValueError:
            _malformed(st, opt, lineno, line, "fields 4-6 are not coordinates")
            continue

        # No element column is not a defect for this format, so it is not
        # reported as one -- but the name is still a better source than the
        # pre-2.1 seven-letter guess.
        element = "" if opt.legacy_elements else infer_element(name, "", True)

        structure.atoms.append(
            Atom(
                x=x, y=y, z=z, name=name, element=element, serial=serial,
                hetatm=True, element_inferred=bool(element),
            )
        )

        if opt.max_atoms is not None and len(structure.atoms) >= opt.max_atoms:
            break

    st.accepted = len(structure.atoms)
    return structure, st


# ----------------------------------------------------------------------
# Types
# ----------------------------------------------------------------------


def assign_types(
    structure: Structure, legacy: bool = False, stats: ParseStats | None = None
) -> None:
    """
    Give every atom the table index that selects its texture.

    Anything with no dedicated texture is :data:`ELEMENT_UNKNOWN` and renders
    as the generic ``Atom_X``, so an unrecognised element cannot vanish from
    a scene the way it did before 2.1.  Under ``legacy`` the pre-2.1 guess
    runs instead and those atoms are dropped by the writer -- reproducing an
    old render means reproducing the omission -- but the drop is counted and
    reported either way.
    """
    st = stats or ParseStats()

    for atom in structure.atoms:
        if legacy:
            symbol = legacy_guess_symbol(atom.name, atom.hetatm)
        else:
            symbol = atom.element

        atom.type_index = element_index(symbol)

        if atom.type_index == ELEMENT_UNKNOWN:
            st.generic += 1
            st.note_generic(atom.element if (atom.element and not legacy) else atom.name)


# ----------------------------------------------------------------------
# Front door
# ----------------------------------------------------------------------


def read_structure(
    path: str, options: ParseOptions | None = None
) -> tuple[Structure, ParseStats]:
    """
    Read a structure from a path, choosing the reader by format.

    ``path`` may be a stem in the historical style (``1crn``), a full path
    with an extension, a compressed file, or ``-`` for standard input.
    Returns the structure with element types already assigned, and the
    statistics describing what the parse skipped or inferred.
    """
    from .mmcif import read_mmcif  # imported here to keep the module optional

    opt = options or ParseOptions()
    resolved = resolve_input_path(path, opt.fmt)

    with open_text(resolved) as handle:
        text = handle.read()

    lines = text.splitlines()
    fmt = opt.fmt
    if fmt is InputFormat.AUTO:
        fmt = detect_format(resolved, lines[:200])

    if fmt is InputFormat.CIF:
        structure, stats = read_mmcif(lines, opt)
    elif fmt is InputFormat.ATM:
        structure, stats = read_atm(lines, opt)
    else:
        structure, stats = read_pdb(lines, opt)

    structure.source = "standard input" if resolved == "-" else resolved
    assign_types(structure, opt.legacy_elements, stats)
    return structure, stats
