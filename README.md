# pdb2pov

**Converts Brookhaven PDB atomic structure files into POV-Ray scenes.**

pdb2pov 2.2 — Copyright (c) 1993-2026 Eric G. Suchanek, Ph.D.
Subject to the GNU License.

Written in 1993 for the Amiga and UNIX. Modernised in 2026: the sources are
prototyped C17, and the scenes it writes are POV-Ray 3.7.

> **There is also a Python port**, `pypdb2pov`, in [`python/`](python/). It
> writes byte-identical scenes from the same command lines, and adds what the
> C cannot do: PDBx/mmCIF and compressed input, model selection, element
> inference that knows more than seven first letters, a linear-time bond
> search, and a library API. Its command is spelled `pypdb2pov`, so it can
> share a `PATH` with the binary this Makefile builds. See
> [`python/README.md`](python/README.md). Everything below describes the C
> program.

---

## Building

```sh
make
```

That is the whole story now. There are no portability flags to arrange, no
force-included prototype headers, and no need to disable `_FORTIFY_SOURCE` —
see [What changed in 2.0](#what-changed-in-20). The build is clean under
`-Wall -Wextra -Wpedantic`.

Other targets:

| Target | Effect |
|--------|--------|
| `make` | build `pdb2pov` |
| `make check` | build, then convert the bundled `1CRN.pdb` three ways |
| `make clean` | remove object files |
| `make distclean` | also remove the binary and `make check` artefacts |

---

## Usage

```
pdb2pov InputFile OutputFile [options]
```

Filenames are given **without extensions** — `.pdb` (or `.atm`) and `.pov`
(or `.inc`) are appended.

```sh
pdb2pov crambin crambin -s -h -b -d 1.5 -x 90
```

Converts `crambin.pdb` to `crambin.pov` with a checkered ground and sky,
rotates the structure 90° about X, and renders ball-and-stick with a 1.5 Å
bond cutoff.

### Options

| Flag | Effect |
|------|--------|
| `-v` | van der Waals radii (default) |
| `-c` | covalent radii |
| `-b` | ball and stick |
| `-q` | ball and stick with glass atoms |
| `-d x.x` | bond cutoff threshold in ångströms (default 2.2) |
| `-r x.x` | scale factor applied to all atomic radii |
| `-o` | object only — write a `.inc` with no camera or lights |
| `-t` | input is in `.atm` format |
| `-p` | plain white sky, no ground |
| `-s` | cloudy sky |
| `-g` | plain ground |
| `-h` | checkered ground |
| `-a` | area light instead of a point light |
| `-x -y -z` | absolute rotation about each axis, in degrees |
| `--chain IDS` | restrict to the given chain IDs, e.g. `--chain AB` |
| `--keep-altlocs` | keep every alternate conformation |
| `--legacy-elements` | guess elements from atom names, pre-2.1 style |

Rendering styles depend on the bundled include files (`atoms2.inc`,
`atoms_vdw.inc`, `atoms_covalent.inc`, `atoms_cpk.inc`, `atoms_glass2.inc`),
which must be findable by POV-Ray — either in the working directory or on its
library path (`+L`).

> `atoms_cpk.inc` ships but is unreachable: no flag selects it, and none did
> in 1993 either. Its radii are identical to `atoms_vdw.inc`'s in any case.

---

## Preparing input

The parser predates several PDB conventions, and earlier releases could read
past the end of short records — which is why trimming a modern file to its
coordinate records first was standard advice. That is no longer required for
correctness: record lengths are checked, and the bundled `1CRN.pdb` converts
to a byte-identical scene whether or not it is trimmed.

Trimming is still worth doing when you want to *choose* what appears — for
instance dropping crystallographic waters:

```sh
grep -E "^(ATOM|END)" 1CRN.pdb > crambin.pdb        # protein only
grep -E "^(ATOM|HETATM|END)" 1CRN.pdb > crambin.pdb # include heteroatoms
```

If a structure misbehaves, stripping column 22 restores the pre-1996 layout
the parser was written against:

```sh
awk '/^ATOM/ {print substr($0,1,21) " " substr($0,23)}' in.pdb > out.pdb
```

### What the parser does with a record

Records are read by fixed PDB column ranges — atom name from 13–16, altLoc
from 17, chain from 22, coordinates from 31–38 / 39–46 / 47–54, element from
77–78. A record that does not honour those columns falls back to the
free-form scan earlier versions used, and the fallback is reported, so
nothing that converted before stops converting.

**Elements come from columns 77–78.** Before 2.1 the element was guessed from
the first one or two characters of the atom name, which is wrong for any
two-letter element sharing a first letter with a one-letter one:

| Record | Element | Pre-2.1 result | 2.1 |
|--------|---------|----------------|-----|
| `NA` | sodium | nitrogen | sodium → `Atom_X` |
| `CL` | chlorine | carbon | chlorine → `Atom_X` |
| `F` | fluorine | iron | fluorine → `Atom_X` |
| `ZN` | zinc | *silently dropped* | zinc → `Atom_X` |
| `CA` (ATOM) | carbon (Cα) | carbon | carbon |
| `CA` (HETATM) | calcium | calcium | calcium |

Files old enough to have no element column fall back to the name guess, and
say so. The guess also required a special case to tell a protein alpha carbon
from calcium; with a real element column that distinction is free, so the
hack is now confined to the legacy path.

**Elements with no dedicated texture render as `Atom_X`,** a neutral grey
sphere, rather than disappearing. Before 2.1 they were dropped without a
message, so an ion could vanish while the header atom count still looked
right. The count and the symbols involved are now reported. See
[Elements](#elements) for what has a texture.

**Alternate conformations are filtered.** Only the blank and `A` altLoc
indicators are kept. Retaining all of them — the pre-2.1 behaviour — gives
overlapping spheres at nearly identical positions plus spurious bonds between
the A and B conformers; a 7-record two-conformer test file yields 7 atoms and
13 bonds where the correct answer is 4 and 3. Use `--keep-altlocs` to keep
them.

**Hydrogen is identified by element, not by name.** The bonding rules cap
hydrogens at one bond; previously any atom whose name began with `H` was
capped, so mercury bridging two cysteines got one bond instead of two.

`--legacy-elements` restores the whole pre-2.1 arrangement — name-based
guessing, the alpha-carbon hack, name-based hydrogen detection, and dropping
unrecognised elements — for reproducing the appearance of an existing render.
The drop is reported rather than silent even then.

Structures made only of H, C, N, O, S and P are unaffected by any of this.
Crambin converts to byte-identical output before and after.

---

## Rendering the output

Scenes are POV-Ray 3.7 and render as generated:

```sh
povray +Icrambin.pov +W800 +H600 +A0.3 -D +L.
```

No `#version` pragma needs prepending, and there are no parse warnings.
Earlier releases wrote POV-Ray 2.x, which POV-Ray 3.5+ rejects outright
unless the language version is pinned; the documented workaround was to
prepend `#version 3.1;` and tolerate the warnings that followed. That is no
longer necessary.

Scenes declare `assumed_gamma 1.0`, the modern default. Pre-3.7 renders were
effectively uncorrected, so colours are lighter than the 1994 originals; add
`global_settings { assumed_gamma 2.2 }` to your host scene if you want the
historical look.

---

## Output format

The header records what was converted, including the bounding geometry:

```
//	Atoms:  327
//	Extent:	Xmin: -14.866 Xmax: 17.515,
//		Ymin: -12.803, Ymax: 13.650
//		Zmin: -15.113 Zmax: 16.889
//	Enclosing Sphere: 18.759
```

The structure is centred on the origin and the camera placed to frame it, so
the enclosing sphere radius and camera distance are all that is needed to
bound the scene's depth — useful when driving the scene through a renderer
that cares about near and far extents.

The enclosing radius is also emitted as a POV-Ray float, so a host scene can
read it directly rather than scraping the comment:

```povray
#declare crambin_enclosing_radius = 18.759;
#declare crambin_obj = union { /* atoms and bonds */ }
#declare crambin     = object { crambin_obj }
```

With `-o` the file is a `.inc` that saves and restores the language version
around its own declarations, so including it will not switch the host scene
to 3.7 behind your back.

---

## Elements

Thirty-three elements have their own radius, colour and texture:

| Group | Elements |
|-------|----------|
| Organic and biological | H, C, N, O, S, P, Se |
| Halogens | F, Cl, Br, I |
| Alkali and alkaline earth | Li, Na, K, Mg, Ca |
| Transition and heavy metals | Mn, Fe, Co, Ni, Cu, Zn, Mo, W, Ag, Cd, Pt, Au, Hg |
| Other | B, Si, As, Xe |

Anything else renders as `Atom_X`, a neutral grey sphere, and is reported.
Deuterium and tritium are drawn as hydrogen.

**Radii.** The eight elements pdb2pov has always drawn keep their 1994 values
(Pauling and the CAChe software). Elements added in 2.2 use Bondi van der
Waals radii and Cordero covalent radii; where a metal's covalent radius
depends on spin state, the low-spin value is used. `atoms_cpk.inc` remains
identical to `atoms_vdw.inc` — Corey-Pauling-Koltun space-filling radii *are*
the van der Waals radii — and is still unreachable, since no flag selects it.

**Colours.** The original eight are not the CPK convention: carbon is green,
phosphorus yellow, iron dark purple, calcium white. Changing them would alter
every existing render, so they stay. Elements added in 2.2 use Jmol/CPK
colours.

The glass textures used by `-q` differ from the solid ones for two elements —
calcium is white solid but dark purple in glass, iron dark purple solid but
gold in glass. That inconsistency dates from the 1994 files and is preserved
for the same reason.

### Adding an element

Four places, and `make check` will tell you if you miss one:

1. a row in the `ELEMENTS` table in `pdb2pov.c`, giving the PDB symbol and the
   POV-Ray identifier suffix;
2. `<SYM>_RAD` in `atoms_vdw.inc`, `atoms_covalent.inc` and `atoms_cpk.inc`;
3. `Atom_<Suffix>` in `atoms2.inc`;
4. `Glass_<Suffix>` and `Atom_Glass_<Suffix>` in `atoms_glass2.inc`.

`make check` converts the bundled `elements.pdb` — one atom of every element
in the table — and renders it. A row added without its declarations is an
undeclared POV-Ray identifier, which is a parse error, so the omission fails
the build rather than surfacing at someone's first metalloprotein.

---

## What changed in 2.2

- Thirty-three elements have dedicated radii, colours and textures, up from
  eight. See [Elements](#elements).
- The element dispatch is a table rather than two parallel switch statements,
  so adding one is a row and four declarations.
- `elements.pdb` and a `make check` step verify the table and the include
  files agree.

Existing scenes are unaffected: the original eight keep their radii and
colours exactly, and crambin converts to byte-identical output.

---

## What changed in 2.1

The parser moved from free-form scanning to the PDB column layout, and three
behaviours that silently produced wrong scenes were corrected. See
[What the parser does with a record](#what-the-parser-does-with-a-record)
above for the detail; in brief:

- Elements are read from columns 77–78 rather than guessed from the atom name.
- Elements with no dedicated texture render as a generic grey `Atom_X`
  instead of being dropped without a message.
- Alternate conformations are filtered to the blank and `A` indicators.
- Hydrogen is identified by element, so mercury is no longer bond-capped.
- `--chain` restricts conversion to named chains.
- The parse reports what it skipped, guessed at, or fell back on.

`--legacy-elements` restores the pre-2.1 behaviour. Only structures
containing elements beyond H, C, N, O, S and P, or carrying alternate
conformations, render differently.

---

## What changed in 2.0

The 1993 sources were K&R C with `#ifdef AMIGA` branches. They only built on
a 64-bit host with a set of Makefile workarounds — a force-included prototype
header, `-std=gnu89`, and `_FORTIFY_SOURCE` disabled. All of that is gone,
because the underlying problems are fixed rather than suppressed.

**Correctness**

- **Pointer truncation.** `pdb2pov.c` called the allocators in `util.c`
  without declaring them. Under K&R rules an undeclared function is assumed
  to return `int`, so every pointer they returned was truncated to 32 bits on
  a 64-bit host — an immediate segfault. Everything is prototyped now, so the
  compiler checks it.
- **Buffer overrun in the date stamp.** The header timestamp was built with
  `sprintf` into a buffer too small to hold it, which is what `_FORTIFY_SOURCE`
  was tripping on. It also printed `tm_year` — years since 1900 — as a
  two-digit field, so 2026 rendered as `126`. Now `strftime` into a checked
  buffer, with a four-digit year.
- **The bond search skipped atom 0.** The inner loop ran `j > 0` where it
  meant `j >= 0`, so the first atom in every file was never bonded to
  anything. Crambin gains exactly one bond, at the N-terminus.
- **Out-of-bounds reads while parsing.** Column offsets up to 70 were indexed
  unconditionally once a record matched `ATOM`/`HETATM`, reading past the
  string terminator on any shorter line. Lengths are checked now.
- **Duplicated bond logic.** Counting and recording bonds were two near-copies
  of the same loop that had already drifted — one used `hypot()`, the other an
  inline `sqrt()`. A single traversal now does both, so they cannot disagree.
- **Invalid POV identifiers.** The declared object name came from the output
  path unfiltered, so `-o out/1crn` emitted `#declare out/1crn` — which does
  not parse. Names are sanitised.
- **Unchecked option arguments.** `-d`, `-x`, `-y`, `-z` and `-r` read past
  the end of `argv` when given as the last word on the line.

**POV-Ray 3.7 output**

- `#declare` of a float now carries the terminating semicolon 3.5+ requires.
- `sphere { <0,0,0> RADIUS }` gained the comma between centre and radius.
- `light_source` is no longer wrapped in `object { }` — in 3.x a light source
  is not an object.
- `refraction` and `ior` moved out of `finish{}`, which no longer accepts
  them, into `interior{}`.
- 2.x `colour_map` ramps (`[start, end colourA colourB]`) became 3.x control
  point lists, and `colour`/`quick_colour` became `color`/`quick_color`.
- `bounded_by { sphere { ... } }` is gone. POV-Ray 3.x bounds CSG
  automatically and warns that the manual sphere is redundant; the radius is
  published as a float instead.
- `right <4/3, 0, 0>` became `right x*image_width/image_height`. At 4:3 these
  are identical; at any other aspect ratio the fixed form rendered with
  non-square pixels.
- Scenes declare `assumed_gamma`, without which POV-Ray 3.7 refuses to render.

**Housekeeping**

- The `#ifdef AMIGA` branches, `asyncio.c` and `m68040.c` are removed; git
  history preserves them.
- `pdb2light.h` and `pdb2pov_protos.h` are consolidated into `pdb2pov.h` and
  `util.h`. The old `pdb2light.h` defined `X`, `Y` and `Z` twice and carried a
  pile of constants nothing referenced.
- Dead code is gone: `flip_yaxis`, `translate_protein`, `compute_max_scale`,
  `inc_rotmat`, `matmul`, `print_coords`, `fgetpdbxyz`, and the unreferenced
  `-m` and `-i` options.
- The allocators return `NULL` on failure instead of calling `exit()` from
  inside the allocator, so callers can release what they already hold. They
  also no longer offset returned pointers by a subscript base that was always
  zero, which was undefined behaviour for any other value.
- The parser no longer echoes every input line to stdout.
- `-lm` is linked explicitly. macOS supplies libm via libSystem, so its
  absence went unnoticed for thirty years; on GNU/Linux the link failed.

**Deliberately unchanged.** Scene *geometry* is identical. The 1993 code
converted degrees to radians with hand-rolled constants (`57.29` in the
rotation matrix, `57.298` for camera and light placement), both truncations
of 180/π. They are kept verbatim because downstream documentation quotes the
figures they produce — crambin's camera distance of 40.075 and enclosing
sphere of 18.759. Correcting them would silently move every camera in every
existing scene, for a difference of well under a tenth of a degree.

---

## See also

- [`python/`](python/) — `pypdb2pov`, the Python port. Same scenes from the
  same flags, plus mmCIF, compressed input, model and conformer selection, a
  linear-time bond search, and an importable API. Its test suite diffs its
  output against this program's, so the two cannot drift apart.
- [quiltwright](https://github.com/suchanek/quiltwright) — renders `.pov`
  scenes as multi-view quilts for Looking Glass holographic displays. Its
  `docs/pdb2pov.md` covers turning structures into holograms, and uses the
  bounding-sphere values above to set the holographic focal plane. Note that
  its build and `#version 3.1;` instructions describe pdb2pov 1.19 and are
  superseded by this release.
- [proteusPy](https://github.com/suchanek/proteusPy) — a modern descendant of
  the same idea: atoms to spheres, bonds to split-coloured cylinders, thirty
  years later and targeting VTK instead of a ray-tracer.
