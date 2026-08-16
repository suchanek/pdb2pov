# pypdb2pov

**Converts Brookhaven PDB and PDBx/mmCIF atomic structure files into POV-Ray
scenes.**

A Python port of [`pdb2pov.c`](../pdb2pov.c) 2.2, in a directory that stands
on its own: no dependencies beyond the standard library, the POV-Ray include
files bundled inside the package, and nothing outside `python/` needed to
build, test or copy it elsewhere.

pypdb2pov 0.1.0 — Copyright (c) 1993-2026 Eric G. Suchanek, Ph.D.
Subject to the GNU License.

> **Two version numbers, deliberately.** `pypdb2pov.__version__` is this
> package's, and starts at 0.1.0 because the package is new.
> `pypdb2pov.PDB2POV_VERSION` is 2.2 — the pdb2pov release whose scenes it
> writes, tracked from `pdb2pov.h`. Scene headers name both:
>
> ```
> // Prepared by pypdb2pov 0.1.0 (pdb2pov 2.2) from 1CRN.pdb on 2026-08-16 …
> ```
>
> The distinction lets a port-only fix ship as 0.1.1 without a scene claiming
> the thirty-year-old C program changed. **The name is also the reason there
> is no collision**: the C program installs a binary called `pdb2pov`, and
> this one installs `pypdb2pov`, so both can be on `PATH` at once.

---

## The promise

**Existing command lines produce the same scenes.** Every flag means what it
has meant since 1993, and the output is byte-identical to the C program's for
the same input — apart from the header's `Prepared by` line, which carries the
timestamp. `tests/test_against_c.py` builds the C program and diffs the two
across nine flag combinations, the 33-element grid, and the legacy element
path, so it stays that way.

That matters because the numbers are quoted downstream: quiltwright's
`docs/pdb2pov.md` pins crambin's camera distance at 40.075 and its enclosing
sphere at 18.759. Both come out of the deliberately truncated
degree-to-radian constants the C uses, and both are reproduced here.

Everything below is what is *new*.

---

## Installing

```sh
cd python
pip install -e .          # or: pip install .
```

Or don't install it at all — the package has no dependencies, so
`PYTHONPATH=python/src python3 -m pypdb2pov ...` works from a checkout, and
copying `python/` somewhere else works too.

Python 3.10 or newer.

---

## Usage

```
pypdb2pov InputFile OutputFile [options]
```

```sh
pypdb2pov crambin crambin -s -h -b -d 1.5 -x 90
```

Filenames may be given without extensions, as they always could — `.pdb` (or
`.cif`, `.ent`, `.atm`, or a compressed variant) is appended to the input and
`.pov`/`.inc` to the output. A path that exists as given is used as given, so
`4hhb.cif.gz` works too. `-` is standard input or standard output.

> **`-h` is the checkered ground, not help.** It has been since 1993, so
> `--help` prints the usage message instead.

### The original flags

| Flag | Effect |
|------|--------|
| `-v` | van der Waals radii (default) |
| `-c` | covalent radii |
| `-b` | ball and stick |
| `-q` | ball and stick with glass atoms |
| `-d x.x` | bond cutoff threshold in ångströms (default 2.2) |
| `-r x.x` | scale factor applied to all atomic radii |
| `-o` | object only — write a `.inc` with no camera or lights |
| `-p` | plain white sky, no ground |
| `-s` | cloudy sky |
| `-g` | plain ground |
| `-h` | checkered ground |
| `-a` | area light instead of a point light |
| `-x -y -z` | absolute rotation about each axis, in degrees |
| `-t` | input is in `.atm` format |
| `--chain IDS` | restrict to the given chain IDs |
| `--keep-altlocs` | keep every alternate conformation |
| `--legacy-elements` | guess elements from atom names, pre-2.1 style |

### What is new

| Flag | Effect |
|------|--------|
| `--format {auto,pdb,cif,atm}` | override format detection |
| `--model N` | convert model N (default: the first in the file) |
| `--all-models` | convert every model at once |
| `--altloc {a,first,occupancy,all}` | how to resolve alternate conformations |
| `--no-hetatm`, `--no-water` | drop heteroatoms or waters without pre-filtering |
| `--bonds {distance,covalent}` | bond by one cutoff, or by covalent radii |
| `--bond-tolerance x.x` | slack added to the radius sum in covalent mode |
| `--strict` | fail on an unparseable record instead of skipping it |
| `--name IDENT` | POV-Ray identifier for the declared objects |
| `--no-timestamp` | omit the date, so two runs are byte-identical |
| `--info` | report what the file contains and write nothing |
| `--include-dir` | print where the bundled `.inc` files live |
| `--list-elements` | print the elements with a dedicated texture |
| `--quiet` | only report problems |

---

## What the port improves

### It reads mmCIF

The PDB format ran out of room years ago — five columns for a serial number,
one for a chain — so anything above 99,999 atoms or 62 chains is distributed
only as PDBx/mmCIF, and `pdb2pov.c` cannot read a word of it. A ribosome, a
capsid or a large complex had no route into a scene at all.

```sh
pypdb2pov 6xyz.cif.gz capsid -b -o
```

The reader handles the `_atom_site` loop with the syntax wwPDB files actually
use: quoted values, semicolon-delimited text blocks, nulls, and multiple
models. `auth_*` columns are preferred over `label_*` where both exist, so
`--chain A` means the chain the literature calls A.

### It opens compressed files

`.gz`, `.bz2` and `.xz` are read transparently, which is how the wwPDB ships
them. No decompression step, no temporary file.

### It understands models

`MODEL`/`ENDMDL` are records, not accidents. The C stopped reading at the
first line beginning `END`, which happens to be `ENDMDL` — the right answer
for an NMR ensemble, reached by the wrong route, with no way to ask for model
7. Here the default is the first model, `--model N` picks another, and
`--all-models` keeps them all.

### It infers elements from more than seven letters

Elements come from columns 77-78 when they are there. When they are not — the
file predates the column — the fallback now reads the atom name the way the
format specifies rather than switching on its first letter:

| Record | Residue | Pre-2.1 guess | Here |
|--------|---------|---------------|------|
| ` NA ` | `HEM` | **sodium** | nitrogen (haem N-alpha) |
| ` CD ` | `GLU` | carbon | carbon (glutamate C-delta) |
| ` CD ` | ligand | carbon | cadmium, if column-aligned as one |
| `ZN  ` | `ZN` | *dropped* | zinc |
| `FE  ` | `HEM` | iron | iron |
| `SE  ` | `MSE` | *dropped* | selenium |

Three rules do the work, in order: an atom name that matches its residue name
and is an element symbol is that element (`ZN`/`ZN`); inside a standard
residue the element is the leading letter, full stop; otherwise the column
alignment decides, because the format right-justifies a two-letter symbol
into columns 13-14.

Where the answer is genuinely ambiguous — a hand-edited ` FE ` sitting in
column 14, where the format says one-letter — the column rule wins and the
conversion **says so**:

```
  1 atom name(s) read by the PDB column rule where a two-letter element was
    also possible ('FE' as F not FE);
      add an element column in 77-78, or use --legacy-elements, to override
```

`--legacy-elements` still reproduces the 1993 behaviour exactly, mistakes
included, for reproducing an existing render.

### It gives alternate conformations a policy

| `--altloc` | Behaviour |
|------------|-----------|
| `a` | keep the blank and `A` indicators — the C's behaviour, and the default |
| `first` | keep whichever indicator comes first for each atom |
| `occupancy` | keep the highest-occupancy conformer, the crystallographer's answer |
| `all` | keep everything, as pdb2pov did before 2.1 |

`a` loses a residue entirely when a depositor labels its conformers `B` and
`C` with no `A`. `first` does not, and is otherwise identical on
conventionally ordered files. In 1CBN that is fourteen atoms — the side
chains of Pro 22 and Leu 25 — that the C drops without a word.

The choice is made **per residue**, which matters for microheterogeneity: in
that same entry residue 22 is modelled as serine at 0.20 occupancy *and*
proline at 0.60, sharing one sequence position. Choosing atom by atom would
take proline's ring and serine's hydroxyl from the same place and draw a
residue that does not exist. `--altloc occupancy` picks the proline, whole.

### The bond search is linear, not quadratic

The C compares every pair of atoms, twice — once to count and once to
record. That is fine for crambin's 327 atoms and hopeless above about ten
thousand. Atoms now go into a uniform grid of cells one cutoff wide, so each
one only looks at the 27 cells around it. 50,000 atoms takes well under a
second.

The results are identical, pair for pair and in the same order: the
traversal order the hydrogen rule depends on is reproduced deliberately, and
`tests/test_bonds.py` diffs the grid search against a transcription of the
C's loop over random structures.

`--bonds covalent` adds what a single cutoff cannot express — a 1.1 Å C-H and
a 2.05 Å disulphide in the same structure, without also bonding every carbon
to its neighbour's neighbour.

### Nothing is dropped silently, and problems have line numbers

```
Scanning atom file <4hhb.pdb>... got <4779> atoms.
  216 atom(s) in alternate conformations skipped (--keep-altlocs to keep them)
  2 unparseable coordinate record(s) skipped
      line 1841: ATOM   1839  CB  SER A 234      ????    ????   ????
  4 atom(s) have no dedicated texture and render as Atom_X (GD, RU)
```

`--strict` turns the skip into a failure. `--info` reports the same census
without writing anything.

### No fixed limits

The C read into a 256-byte line buffer and a pre-counted fixed-size array; a
line longer than the buffer was silently split into two records. There are no
buffers here, and no pre-count pass — files are read once.

### It is a library

```python
from pypdb2pov import convert

convert("1crn.pdb", "crambin", ball_stick=True, bond_threshold=1.9)
```

or, with the pieces exposed:

```python
from pypdb2pov import (
    ParseOptions, SceneOptions, read_structure,
    find_bonds, prepare_structure, write_scene,
)

structure, stats = read_structure("4hhb.cif.gz", ParseOptions(chains="A"))
print(structure.title, len(structure), structure.element_counts())

options = SceneOptions(ball_stick=True, object_only=True, name="hemoglobin_a")
prepare_structure(structure, options)                    # rotate, centre, flip
bonds = find_bonds(structure, options.bond_threshold)
write_scene(structure, options, "hemoglobin_a.inc", bonds)
```

`Structure` is a plain list of `Atom` dataclasses carrying everything the
readers recover — name, element, residue, chain, insertion code, altLoc,
occupancy, B-factor, formal charge, model — so it is usable as a PDB/mmCIF
reader in its own right, not only as a step towards a scene.

---

## Rendering

Scenes reference the bundled include files by name, so POV-Ray needs the
directory on its library path:

```sh
povray +Icrambin.pov +W800 +H600 +A0.3 -D +L"$(pypdb2pov --include-dir)"
```

The includes are shipped inside the package, which is what makes this
directory self-contained. `tests/test_against_c.py` asserts they are still
byte-identical to the copies in the C tree, so the two cannot drift.

---

## Testing

```sh
make test          # the whole suite
make check         # the suite, plus a reminder of what was skipped
```

or directly:

```sh
python -m pytest
```

112 tests. Those that compare against the C program build it first and skip
cleanly when the C sources are not alongside the package — so a copy of
`python/` on its own still tests green.

---

## Differences from the C, in full

Everything here is deliberate. Nothing else differs.

| Area | C 2.2 | Python |
|------|-------|--------|
| Scene output | — | byte-identical, except the `Prepared by` line, which also names this package and the source file |
| Element inference with no element column | seven-letter first-character guess | column- and residue-aware, with ambiguity reported |
| `.atm` element handling | first-character guess | the same inference as PDB; `--legacy-elements` restores the guess |
| Formats | PDB, `.atm` | plus mmCIF, plus `.gz`/`.bz2`/`.xz`, plus stdin |
| Models | first `END`-prefixed line ends the read | `MODEL`/`ENDMDL`, `--model`, `--all-models` |
| Bond search | O(N²), twice | O(N) via a grid, same results and order |
| Line length | 256 bytes, longer lines split | unbounded |
| Atom count | pre-counted, fixed array | unbounded |
| Exit statuses | 0/2/3/4/5/6 | the same |
| Command name | `pdb2pov` | `pypdb2pov`, so both can share a `PATH` |

The `-t` `.atm` path is the one place a *default* differs in output: the
format has no element column, and the port infers elements there the same way
it does for PDB rather than falling back to the 1993 guess. Pass
`--legacy-elements` for the old behaviour.

---

## See also

- [`../README.md`](../README.md) — the C program, the element table, and the
  history of what changed in 2.0, 2.1 and 2.2.
- [quiltwright](https://github.com/suchanek/quiltwright) — renders these
  scenes as multi-view quilts for Looking Glass holographic displays.
- [proteusPy](https://github.com/suchanek/proteusPy) — the same idea again,
  targeting VTK.
