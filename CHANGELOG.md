# Changelog

All notable changes to `pdb2pov` are documented here.

The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/).
Versions are the C program's own, as reported by `PDB2POV_VERSION` in
`pdb2pov.h`, and predate this file: the entries below 2.2 were reconstructed
from the README sections they replace and from the RCS logs, which are dated
1993-94.

## [Unreleased]

### Removed

- **The `python/` tree.** `pypdb2pov` now lives in its own repository at
  <https://github.com/Flux-Frontiers/pypdb2pov> and is released from there.
  It was carried here while it was a port of this program rather than a
  package in its own right; keeping a second copy alongside the C would only
  invite the two to drift.

  This repository is again what it was: the C program, its include files, and
  its documentation. Nothing about the C changes, and the scenes both write
  are still byte-identical.

  **One thing does not survive the move.** `python/tests/test_against_c.py`
  built this program and diffed its output against the port's, across the
  flag matrix, the element grid and the bundled include files. That suite
  needed the two to sit in one tree, which is exactly what has ended, and the
  new repository does not carry it. The byte-identical claim is therefore no
  longer enforced by anything on either side. Re-establishing it means a
  conftest in the new repository that builds this program when it can find
  it and skips when it cannot.

## [2.2] - 2026

### Added

- Thirty-three elements have dedicated radii, colours and textures, up from
  eight.
- `elements.pdb` and a `make check` step verify that the element table and
  the include files agree.

### Changed

- Element dispatch is a table rather than two parallel `switch` statements,
  so adding an element is a row and four declarations.

Existing scenes are unaffected: the original eight keep their radii and
colours exactly, and crambin converts to byte-identical output.

## [2.1] - 2026

### Changed

- The parser reads the PDB **column layout** rather than scanning free-form.
  Elements come from columns 77-78 instead of being guessed from the atom
  name, which was wrong for any two-letter element sharing a first letter
  with a one-letter one.
- Hydrogen is identified by element, so mercury is no longer bond-capped.
- Alternate conformations are filtered to the blank and `A` indicators.
- The parse reports what it skipped, guessed at, or fell back on.

### Added

- `--chain` restricts conversion to named chains.
- `--legacy-elements` restores the pre-2.1 name guessing, for reproducing an
  old render exactly.

### Fixed

- An element with no dedicated texture renders as a generic grey `Atom_X`
  instead of being dropped without a message.

## [2.0] - 2026

### Fixed

- **Pointer truncation.** The allocators in `util.c` were called without
  being declared, so under K&R rules every pointer they returned was
  truncated to 32 bits on a 64-bit host -- an immediate segfault. Everything
  is prototyped now.
- **Buffer overrun in the date stamp.** The header timestamp was built with
  `sprintf` into a buffer too small for it, and printed `tm_year` as a
  two-digit field, so 2026 rendered as `126`. Now `strftime` into a checked
  buffer, with a four-digit year.
- **The bond search skipped atom 0.** The inner loop ran `j > 0` where it
  meant `j >= 0`, so the first atom in every file was never bonded. Crambin
  gains exactly one bond, at the N-terminus.
- **Out-of-bounds reads while parsing.** Column offsets up to 70 were indexed
  unconditionally once a record matched `ATOM`/`HETATM`, reading past the
  terminator on any shorter line. Lengths are checked now.

### Changed

- The 1993 K&R sources with `#ifdef AMIGA` branches are now prototyped C17,
  and the scenes written are POV-Ray 3.7.
- The build needs no workarounds: no force-included prototype header, no
  `-std=gnu89`, no disabling `_FORTIFY_SOURCE`. It is clean under `-Wall
  -Wextra -Wpedantic`.

## [1.19] - 1994

The last release of the original line, written for the Amiga and UNIX. Its
scenes are POV-Ray 1.x/2.x and need `+MV3.1` or a `#version` prepend to parse
under 3.5 and later.

[Unreleased]: https://github.com/suchanek/pdb2pov/compare/v2.2...HEAD
