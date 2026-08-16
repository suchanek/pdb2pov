"""
A PDBx/mmCIF reader, sufficient for coordinates.

This is the reason the Python port exists as more than a translation.  The
PDB format ran out of room years ago -- five columns for a serial number,
one for a chain -- so anything above 99,999 atoms or 62 chains is
distributed only as mmCIF, and ``pdb2pov.c`` cannot read a word of it.  A
ribosome, a capsid or a large complex simply had no route into a scene.

Only what a picture needs is parsed: the ``_atom_site`` loop, plus the entry
title if one is present.  Reading the whole dictionary would mean a CIF
object model, and nothing downstream would look at it.

The tokeniser handles the parts of the syntax that appear in wwPDB files:
quoted values, semicolon-delimited text blocks, ``loop_`` tables, key/value
pairs, and ``.``/``?`` for null.  Multi-line values inside a loop are rare
but legal, and are handled.
"""

from __future__ import annotations

from typing import Iterable, Iterator

from .elements import infer_element_verbose, normalise_symbol
from .options import ALL_MODELS, AltLocPolicy, ParseOptions
from .readers import apply_altloc_policy  # shared with the PDB reader
from .structure import WATER_RESIDUES, Atom, ParseStats, Structure

__all__ = ["read_mmcif", "tokenise"]

_NULL = (".", "?")

#: Unquoted tokens that are syntax, never data.  A key whose value is missing
#: must not swallow the next one -- a truncated ``_struct.title`` followed by
#: ``loop_`` would otherwise eat the keyword and hide the atoms behind it.
_RESERVED = ("loop_", "stop_", "global_")


def _is_syntax(token: str, quoted: bool) -> bool:
    if quoted:
        return False
    lowered = token.lower()
    return (
        token.startswith("_")
        or lowered in _RESERVED
        or lowered.startswith("data_")
        or lowered.startswith("save_")
    )


def _clean(value: str) -> str:
    """CIF nulls become the empty string."""
    return "" if value in _NULL else value


def tokenise(lines: Iterable[str]) -> Iterator[tuple[str, bool]]:
    """
    Yield ``(token, was_quoted)`` for a CIF stream.

    ``was_quoted`` matters because a quoted ``loop_`` is data, not the
    keyword, and an unquoted ``.`` is null while a quoted one is a full stop.
    """
    pending_text: list[str] | None = None

    for raw in lines:
        line = raw.rstrip("\n").rstrip("\r")

        # Semicolon-delimited text block: a ';' in column 1 opens and closes.
        if pending_text is not None:
            if line.startswith(";"):
                yield ("\n".join(pending_text), True)
                pending_text = None
            else:
                pending_text.append(line)
            continue

        if line.startswith(";"):
            pending_text = [line[1:]]
            continue

        stripped = line.strip()
        if not stripped or stripped.startswith("#"):
            continue

        index = 0
        length = len(line)
        while index < length:
            char = line[index]
            if char in " \t":
                index += 1
                continue
            if char == "#":
                break  # comment to end of line
            if char in "'\"":
                quote = char
                index += 1
                start = index
                # A quote only closes when followed by whitespace or the end
                # of the line, which is what lets O5' appear unquoted inside
                # a quoted value.
                while index < length:
                    if line[index] == quote and (
                        index + 1 >= length or line[index + 1] in " \t"
                    ):
                        break
                    index += 1
                yield (line[start:index], True)
                index += 1
                continue

            start = index
            while index < length and line[index] not in " \t":
                index += 1
            yield (line[start:index], False)

    if pending_text is not None:  # unterminated block; take what there is
        yield ("\n".join(pending_text), True)


#: The ``_atom_site`` columns worth reading, in preference order.  The
#: ``auth_*`` variants are the ones that match the PDB-format file for the
#: same entry, so they are preferred where both exist -- a user asking for
#: ``--chain A`` means the chain the literature calls A.
_FIELD_ALIASES = {
    "group": ("group_pdb",),
    "serial": ("id",),
    "name": ("auth_atom_id", "label_atom_id"),
    "alt_loc": ("label_alt_id", "auth_alt_id"),
    "res_name": ("auth_comp_id", "label_comp_id"),
    "chain": ("auth_asym_id", "label_asym_id"),
    "res_seq": ("auth_seq_id", "label_seq_id"),
    "i_code": ("pdbx_pdb_ins_code",),
    "x": ("cartn_x",),
    "y": ("cartn_y",),
    "z": ("cartn_z",),
    "occupancy": ("occupancy",),
    "b": ("b_iso_or_equiv",),
    "element": ("type_symbol",),
    "charge": ("pdbx_formal_charge",),
    "model": ("pdbx_pdb_model_num",),
}


def read_mmcif(
    lines: Iterable[str],
    options: ParseOptions | None = None,
    stats: ParseStats | None = None,
) -> tuple[Structure, ParseStats]:
    """
    Read the first ``_atom_site`` loop of an mmCIF stream.

    Returns the structure and the parse statistics, in the same shape as
    :func:`pypdb2pov.readers.read_pdb`, so callers do not care which format
    they were handed.
    """
    opt = options or ParseOptions()
    st = stats or ParseStats()
    structure = Structure(source_format="cif")

    tokens = tokenise(lines)
    title = ""
    entry_id = ""
    want_all_models = opt.model == ALL_MODELS

    # One token of pushback: a key with no value has to be able to hand the
    # keyword it found back to the main loop instead of consuming it.
    pushback: tuple[str, bool] | None = None

    while True:
        if pushback is not None:
            token, pushback = pushback, None
        else:
            token = next(tokens, None)
        if token is None:
            break
        text, quoted = token

        if not quoted and text.lower() == "loop_":
            headers: list[str] = []
            nxt: tuple[str, bool] | None = None
            while True:
                nxt = next(tokens, None)
                if nxt is None:
                    break
                value, was_quoted = nxt
                if not was_quoted and value.startswith("_"):
                    headers.append(value.lower())
                else:
                    break

            if not headers or not headers[0].startswith("_atom_site."):
                continue

            # `nxt` holds the first data token of the loop.
            _read_atom_site(
                structure, st, opt, tokens, headers, _column_map(headers), nxt,
                want_all_models,
            )
            continue

        # _struct.title is the descriptive one; _entry.id is the four-letter
        # code, and only worth using when there is no title.
        if not quoted and text.lower() in ("_struct.title", "_entry.id"):
            nxt = next(tokens, None)
            if nxt is None:
                continue
            if _is_syntax(*nxt):
                pushback = nxt  # the key had no value; do not eat a keyword
                continue
            value = _clean(nxt[0]).strip()
            if text.lower() == "_struct.title":
                title = title or value
            else:
                entry_id = entry_id or value

    structure.title = title or entry_id
    apply_altloc_policy(structure, opt, st)
    st.accepted = len(structure.atoms)
    return structure, st


def _column_map(headers: list[str]) -> dict[str, int]:
    """Map our field names onto positions in this loop's header list."""
    short = [h.split(".", 1)[1] if "." in h else h for h in headers]
    positions = {name: i for i, name in enumerate(short)}

    columns: dict[str, int] = {}
    for field, aliases in _FIELD_ALIASES.items():
        for alias in aliases:
            if alias in positions:
                columns[field] = positions[alias]
                break
    return columns


def _read_atom_site(
    structure: Structure,
    st: ParseStats,
    opt: ParseOptions,
    tokens: Iterator[tuple[str, bool]],
    headers: list[str],
    columns: dict[str, int],
    first_token: tuple[str, bool] | None,
    want_all_models: bool,
) -> None:
    """Consume the rows of an ``_atom_site`` loop and append the atoms."""
    width = len(headers)
    if width == 0 or "x" not in columns:
        return

    row: list[str] = []
    ordinal = 0
    first_model: int | None = None
    pending: tuple[str, bool] | None = first_token

    while True:
        if pending is not None:
            token, quoted = pending
            pending = None
        else:
            nxt = next(tokens, None)
            if nxt is None:
                break
            token, quoted = nxt

        # A bare keyword or a new data block ends the loop.
        if _is_syntax(token, quoted):
            break

        row.append(token)
        if len(row) < width:
            continue

        ordinal += 1
        if first_model is None:
            # Fixed from the first data row, not the first *accepted* one:
            # a chain filter that rejects row 1 must not leave every model
            # looking like the wanted one.
            first_model = _row_model(columns, row)
        _append_row(structure, st, opt, columns, row, ordinal, first_model, want_all_models)
        row = []

        if opt.max_atoms is not None and len(structure.atoms) >= opt.max_atoms:
            return

    if row:
        st.note_malformed(ordinal + 1, " ".join(row)[:72])


def _row_model(columns: dict[str, int], row: list[str]) -> int:
    index = columns.get("model")
    if index is None:
        return 1
    try:
        return int(_clean(row[index]) or 1)
    except (ValueError, IndexError):
        return 1


def _append_row(
    structure: Structure,
    st: ParseStats,
    opt: ParseOptions,
    columns: dict[str, int],
    row: list[str],
    ordinal: int,
    first_model: int | None,
    want_all_models: bool,
) -> None:
    def get(field: str, default: str = "") -> str:
        index = columns.get(field)
        if index is None:
            return default
        return _clean(row[index])

    group = get("group", "ATOM").upper()
    hetatm = group.startswith("HETATM")

    model_text = get("model")
    try:
        model = int(model_text) if model_text else 1
    except ValueError:
        model = 1
    if model not in st.models_seen:
        st.models_seen.append(model)

    wanted = opt.model if opt.model is not None else (
        first_model if first_model is not None else model
    )
    if not want_all_models and model != wanted:
        st.skipped_model += 1
        return

    alt_loc = get("alt_loc") or " "
    if opt.altloc is AltLocPolicy.A and alt_loc not in (" ", "A", "a"):
        st.skipped_altloc += 1
        return

    chain = get("chain") or " "
    if not opt.accepts_chain(chain):
        st.skipped_chain += 1
        return

    if hetatm and not opt.keep_hetatm:
        st.skipped_hetatm += 1
        return

    res_name = get("res_name")
    if not opt.keep_water and res_name.upper() in WATER_RESIDUES:
        st.skipped_water += 1
        return

    try:
        x = float(get("x"))
        y = float(get("y"))
        z = float(get("z"))
    except ValueError:
        if opt.strict:
            raise
        st.note_malformed(ordinal, " ".join(row)[:72])
        return

    name = get("name")
    element = normalise_symbol(get("element"))
    inferred = False
    if not element:
        st.no_element_column += 1
        # mmCIF atom names are not column-aligned, so the alignment rule in
        # infer_element() cannot fire; the residue-name rules still can.
        element, rejected = infer_element_verbose(name, res_name, hetatm)
        if element:
            inferred = True
            st.inferred_elements += 1
        if rejected:
            st.note_ambiguous(name, element, rejected)

    def number(field: str, default: float) -> float:
        text = get(field)
        if not text:
            return default
        try:
            return float(text)
        except ValueError:
            return default

    try:
        res_seq = int(get("res_seq") or 0)
    except ValueError:
        res_seq = 0

    try:
        serial = int(get("serial") or 0)
    except ValueError:
        serial = ordinal

    structure.atoms.append(
        Atom(
            x=x,
            y=y,
            z=z,
            name=name,
            element=element,
            res_name=res_name,
            chain_id=chain,
            res_seq=res_seq,
            i_code=get("i_code") or " ",
            alt_loc=alt_loc,
            serial=serial,
            occupancy=number("occupancy", 1.0),
            temp_factor=number("b", 0.0),
            charge=number("charge", 0.0),
            hetatm=hetatm,
            model=model,
            element_inferred=inferred,
        )
    )
