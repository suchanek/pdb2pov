# pdb2pov

**Converts Brookhaven PDB atomic structure files into POV-Ray scenes.**

pdb2pov V1.19 — Copyright (c) 1993-2026 Eric G. Suchanek, Ph.D.
Subject to the GNU License.

Written in 1993 for the Amiga and UNIX. It still builds and runs on modern
64-bit hosts, unmodified, given the portability flags below.

---

## Building

```sh
make pdb2pov
```

The Makefile supplies a `PORTFLAGS` set that lets the original sources compile
under a current toolchain without editing them. Three things need handling on
a modern machine:

**Pointer truncation (the one that actually breaks it).** `pdb2pov.c` calls the
allocators defined in `util.c` — `dmatrix`, `dvector`, `cmatrix`, `ivector`,
`imatrix` — without declaring them. Under K&R rules an undeclared function is
assumed to return `int`, so on a 64-bit host every pointer they hand back is
truncated to 32 bits and the program segfaults immediately. `pdb2pov_protos.h`
supplies the missing prototypes and is force-included from the Makefile. The
non-Amiga include path also omits `stdlib.h`, which truncates `malloc` the same
way.

**Fortified libc.** Modern `_FORTIFY_SOURCE` traps some `sprintf` overruns that
were benign in 1993, so it is disabled along with the stack protector.

**C dialect.** `-std=gnu89` for the K&R function definitions, plus
`-Wno-implicit-function-declaration` for several in-file helpers used before
declaration. Those all return `int`, so their implicit declarations are
genuinely harmless.

Both `-O` and `-O0` produce a working binary.

> Object files and a prebuilt binary used to be committed here. They have been
> removed: `make` treated the stale `.o` files as current and linked them
> instead of recompiling, so a fresh clone silently produced a binary of the
> wrong architecture.

---

## Usage

```
pdb2pov InputFile OutputFile [options]
```

Filenames are given **without extensions** — `.pdb` and `.pov` are appended.

```sh
pdb2pov crambin crambin -s -h -b -d 1.5 -x 90
```

Converts `crambin.pdb` to `crambin.pov` with a checkered ground and sky,
rotates the structure 90° about X, and renders ball-and-stick with a 1.5 Å bond
cutoff.

### Options

| Flag | Effect |
|------|--------|
| `-v` | van der Waals radii |
| `-c` | covalent radii |
| `-b` | ball and stick |
| `-q` | ball and stick with glass atoms |
| `-d x.x` | bond cutoff threshold in ångströms |
| `-o` | object only — no camera or lights, for composing into another scene |
| `-t` | input is in `.atm` format |
| `-p` | no sky or ground |
| `-s` | cloudy sky |
| `-g` | plain ground |
| `-h` | checkered ground |
| `-a` | area light |
| `-x -y -z` | absolute rotation about each axis, in degrees |

Rendering styles depend on the bundled include files (`atoms2.inc`,
`atoms_cpk.inc`, `atoms_vdw.inc`, `atoms_covalent.inc`, `atoms_glass2.inc`),
which must be findable by POV-Ray — either in the working directory or on its
library path (`+L`).

---

## Preparing input

The parser predates several PDB conventions. Modern files generally need
trimming to coordinate records:

```sh
grep -E "^(ATOM|HETATM|END)" 1CRN.pdb > crambin.pdb
```

`REMARK 290` crystallographic records in particular will crash it. If a
structure still misbehaves, stripping column 22 restores the pre-1996 layout
the parser was written against:

```sh
awk '/^ATOM/ {print substr($0,1,21) " " substr($0,23)}' in.pdb > out.pdb
```

---

## Rendering the output

The generated scenes use POV-Ray 2.x syntax, where `#declare` statements carry
no trailing semicolon. POV-Ray 3.5 and later reject that unless the language
version is pinned, so prepend one line:

```sh
printf '#version 3.1;\n' | cat - crambin.pov > scene.pov
povray +Iscene.pov +W800 +H600 +FN +A0.3 -D
```

POV-Ray still warns about the missing semicolons. The render is correct.

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

---

## See also

- [quiltwright](https://github.com/suchanek/quiltwright) — renders `.pov`
  scenes as multi-view quilts for Looking Glass holographic displays. Its
  `docs/pdb2pov.md` covers turning structures into holograms, and uses the
  bounding-sphere values above to set the holographic focal plane.
- [proteusPy](https://github.com/suchanek/proteusPy) — a modern descendant of
  the same idea: atoms to spheres, bonds to split-coloured cylinders, thirty
  years later and targeting VTK instead of a ray-tracer.
