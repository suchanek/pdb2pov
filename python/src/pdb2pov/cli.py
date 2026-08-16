"""
Command line interface.

The 1993 grammar is preserved exactly -- ``pdb2pov InputFile OutputFile
[options]``, filenames without extensions, and every single-letter flag
meaning what it has always meant.  In particular **-h is the checkered
ground, not help**, which is why the parser is built with ``add_help=False``
and offers ``--help`` instead.  Anything new is spelled as a long option, so
an existing command line cannot change meaning.
"""

from __future__ import annotations

import argparse
import sys
from typing import Sequence

from . import __version__, include_dir
from .bonds import find_bonds
from .elements import ELEMENTS
from .options import (
    ALL_MODELS,
    AltLocPolicy,
    Backdrop,
    BondMode,
    Ground,
    InputFormat,
    ParseOptions,
    RadiiSet,
    SceneOptions,
)
from .readers import ParseError, read_structure, resolve_input_path
from .scene import PDB2POV_VERSION, POV_VERSION, pov_identifier, prepare_structure, write_scene

__all__ = ["main", "build_parser"]

# Exit statuses, matching the C so scripts wrapping either program agree.
EXIT_OK = 0
EXIT_PARSE_ARGS = 2
EXIT_NO_ATOMS = 3
EXIT_NO_BONDS = 4
EXIT_CANT_WRITE = 5
EXIT_CANT_READ = 6

_EPILOG = """\
examples:
  pdb2pov crambin crambin -s -h -b -d 1.5 -x 90
      crambin.pdb -> crambin.pov, checkered ground and cloudy sky, rotated
      90 degrees about X, ball and stick with a 1.5 angstrom bond cutoff.

  pdb2pov 4hhb.cif.gz hemoglobin -b -o --chain A
      one subunit of a compressed mmCIF entry as a camera-less include.

  pdb2pov 1crn - -b -o > scene.inc
      write the scene to standard output.

note:
  -h is the checkered ground, as it has been since 1993.  Use --help for
  this message.
"""


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="pdb2pov",
        add_help=False,  # -h is the checkered ground
        allow_abbrev=False,
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=(
            f"pdb2pov {PDB2POV_VERSION} -- PDB and mmCIF to POV-Ray "
            f"{POV_VERSION} scene conversion."
        ),
        epilog=_EPILOG,
        usage="%(prog)s InputFile OutputFile [options]",
    )

    parser.add_argument(
        "input",
        nargs="?",
        metavar="InputFile",
        help="input structure; an extension is appended if the path has none, "
        "and '-' reads standard input",
    )
    parser.add_argument(
        "output",
        nargs="?",
        metavar="OutputFile",
        help="output scene; '.pov' or '.inc' is appended, and '-' writes to "
        "standard output",
    )

    rep = parser.add_argument_group("representation")
    rep.add_argument("-v", "-V", dest="vdw", action="store_true",
                     help="van der Waals radii (default)")
    rep.add_argument("-c", "-C", dest="covalent", action="store_true",
                     help="covalent radii")
    rep.add_argument("-b", "-B", dest="ball_stick", action="store_true",
                     help="ball and stick")
    rep.add_argument("-q", "-Q", dest="glass", action="store_true",
                     help="ball and stick with glass atoms")
    rep.add_argument("-d", "-D", dest="bond_threshold", type=float, metavar="X.X",
                     help="bond cutoff threshold in angstroms (default 2.2)")
    rep.add_argument("-r", "-R", dest="radii_scale", type=float, metavar="X.X",
                     help="scale factor applied to all atomic radii")

    scene = parser.add_argument_group("scene")
    scene.add_argument("-o", "-O", dest="object_only", action="store_true",
                       help="object only -- write a .inc with no camera or lights")
    scene.add_argument("-p", "-P", dest="plain_sky", action="store_true",
                       help="plain white sky, no ground")
    scene.add_argument("-s", "-S", dest="sky", action="store_true", help="cloudy sky")
    scene.add_argument("-g", "-G", dest="ground", action="store_true", help="plain ground")
    scene.add_argument("-h", "-H", dest="checker", action="store_true",
                       help="checkered ground (not help; see --help)")
    scene.add_argument("-a", "-A", dest="area_light", action="store_true",
                       help="area light instead of a point light")

    orient = parser.add_argument_group("orientation")
    orient.add_argument("-x", "-X", dest="xrot", type=float, default=0.0, metavar="DEG",
                        help="absolute rotation about X")
    orient.add_argument("-y", "-Y", dest="yrot", type=float, default=0.0, metavar="DEG",
                        help="absolute rotation about Y")
    orient.add_argument("-z", "-Z", dest="zrot", type=float, default=0.0, metavar="DEG",
                        help="absolute rotation about Z")

    inp = parser.add_argument_group("input")
    inp.add_argument("-t", "-T", dest="atm", action="store_true",
                     help="input is in .atm format rather than PDB")
    inp.add_argument("--format", dest="fmt", choices=[f.value for f in InputFormat],
                     default=InputFormat.AUTO.value,
                     help="override format detection (default: auto)")
    inp.add_argument("--chain", dest="chains", metavar="IDS",
                     help="restrict to the given chain IDs, e.g. --chain AB, or a "
                          "comma-separated list for mmCIF's longer names")
    inp.add_argument("--model", dest="model", type=int, metavar="N",
                     help="convert model N (default: the first model in the file)")
    inp.add_argument("--all-models", dest="all_models", action="store_true",
                     help="convert every model at once")
    inp.add_argument("--altloc", dest="altloc",
                     choices=[p.value for p in AltLocPolicy], default=None,
                     help="alternate conformation policy: a (blank or A, the "
                          "default), first, occupancy, all")
    inp.add_argument("--keep-altlocs", dest="keep_altlocs", action="store_true",
                     help="keep every alternate conformation; same as --altloc all")
    inp.add_argument("--no-hetatm", dest="no_hetatm", action="store_true",
                     help="skip HETATM records")
    inp.add_argument("--no-water", dest="no_water", action="store_true",
                     help="skip waters")
    inp.add_argument("--legacy-elements", dest="legacy_elements", action="store_true",
                     help="guess elements from atom names in the pre-2.1 way, and "
                          "drop what the guess cannot place")
    inp.add_argument("--strict", dest="strict", action="store_true",
                     help="fail on an unparseable coordinate record instead of "
                          "skipping it")

    bond = parser.add_argument_group("bonds")
    bond.add_argument("--bonds", dest="bond_mode", choices=[m.value for m in BondMode],
                      default=BondMode.DISTANCE.value,
                      help="distance: one cutoff for every pair (default); "
                           "covalent: the sum of the two covalent radii")
    bond.add_argument("--bond-tolerance", dest="bond_tolerance", type=float,
                      default=0.45, metavar="X.X",
                      help="slack added to the radius sum in covalent mode "
                           "(default 0.45)")

    misc = parser.add_argument_group("other")
    misc.add_argument("--name", dest="name", metavar="IDENT",
                      help="POV-Ray identifier for the declared objects "
                           "(default: derived from OutputFile)")
    misc.add_argument("--no-timestamp", dest="no_timestamp", action="store_true",
                      help="leave the date out of the header, so two runs of the "
                           "same conversion are byte-identical")
    misc.add_argument("--info", dest="info", action="store_true",
                      help="report what the file contains and write nothing")
    misc.add_argument("--include-dir", dest="show_include_dir", action="store_true",
                      help="print the directory holding the bundled .inc files")
    misc.add_argument("--list-elements", dest="list_elements", action="store_true",
                      help="print the elements with a dedicated texture")
    # No short form: -q is the glass representation and has been since 1993.
    misc.add_argument("--quiet", dest="quiet", action="store_true",
                      help="only report problems")
    misc.add_argument("--version", action="version",
                      version=f"pdb2pov {PDB2POV_VERSION} (python {__version__})")
    misc.add_argument("--help", action="help", help="show this message and exit")

    return parser


def _parse_options_from(args: argparse.Namespace) -> ParseOptions:
    fmt = InputFormat(args.fmt)
    if args.atm and fmt is InputFormat.AUTO:
        fmt = InputFormat.ATM

    if args.altloc is not None:
        altloc = AltLocPolicy(args.altloc)
    elif args.keep_altlocs:
        altloc = AltLocPolicy.ALL
    else:
        altloc = AltLocPolicy.A

    chains: str | list[str] | None = None
    if args.chains:
        chains = args.chains.split(",") if "," in args.chains else args.chains

    model = ALL_MODELS if args.all_models else args.model

    return ParseOptions(
        fmt=fmt,
        chains=chains,
        model=model,
        altloc=altloc,
        keep_hetatm=not args.no_hetatm,
        keep_water=not args.no_water,
        legacy_elements=args.legacy_elements,
        strict=args.strict,
    )


def _scene_options_from(args: argparse.Namespace, ident: str) -> SceneOptions:
    radii = RadiiSet.COVALENT if args.covalent else RadiiSet.VDW

    backdrop = Backdrop.NONE
    ground = Ground.NONE
    # -p sets a plain sky and clears the ground, as in the C, where the two
    # assignments happen in the same case.  A later -s or -g overrides.
    if args.plain_sky:
        backdrop = Backdrop.PLAIN
    if args.sky:
        backdrop = Backdrop.SKY
    if args.ground:
        ground = Ground.PLAIN
    if args.checker:
        ground = Ground.CHECKER

    radii_scale = 1.0
    if args.radii_scale is not None:
        radii_scale = args.radii_scale if args.radii_scale > 0.0 else 1.0

    bond_mode = BondMode(args.bond_mode)
    if args.bond_threshold is not None:
        threshold: float | None = args.bond_threshold
    elif bond_mode is BondMode.COVALENT:
        # Let the radii decide rather than clamping them to the historical
        # 2.2 A cutoff, which would cut off exactly the long metal-ligand
        # bonds covalent mode exists to get right.
        threshold = None
    else:
        threshold = 2.2

    return SceneOptions(
        radii=radii,
        backdrop=backdrop,
        ground=ground,
        ball_stick=args.ball_stick or args.glass,
        glass_atoms=args.glass,
        area_light=args.area_light,
        object_only=args.object_only,
        radii_scale=radii_scale,
        bond_threshold=threshold,
        bond_mode=bond_mode,
        bond_tolerance=args.bond_tolerance,
        xrot=args.xrot,
        yrot=args.yrot,
        zrot=args.zrot,
        name=ident,
        timestamp=not args.no_timestamp,
        legacy_drop_unknown=args.legacy_elements,
    )


def main(argv: Sequence[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)

    if args.show_include_dir:
        print(include_dir())
        return EXIT_OK

    if args.list_elements:
        for element in ELEMENTS:
            print(f"{element.symbol:<3} Atom_{element.pov_suffix}")
        print(f"{len(ELEMENTS)} elements with a dedicated texture; "
              "anything else renders as Atom_X")
        return EXIT_OK

    if not args.input:
        parser.print_usage(sys.stderr)
        print("pdb2pov: an input file is required", file=sys.stderr)
        return EXIT_PARSE_ARGS

    if not args.output and not args.info:
        parser.print_usage(sys.stderr)
        print("pdb2pov: an output file is required (or --info)", file=sys.stderr)
        return EXIT_PARSE_ARGS

    popt = _parse_options_from(args)

    # Progress goes to stderr when the scene itself is on stdout.
    stream = sys.stderr if args.output == "-" else sys.stdout

    def say(message: str, end: str = "\n") -> None:
        if not args.quiet:
            print(message, end=end, file=stream, flush=True)

    resolved = resolve_input_path(args.input, popt.fmt)
    say(f"Scanning atom file <{resolved}>... ", end="")

    try:
        structure, stats = read_structure(args.input, popt)
    except ParseError as exc:
        print(f"pdb2pov: {exc}", file=sys.stderr)
        return EXIT_NO_ATOMS
    except OSError as exc:
        print(f"pdb2pov: can't read '{resolved}': {exc.strerror or exc}", file=sys.stderr)
        return EXIT_CANT_READ

    if not structure.atoms:
        say("")
        print(_no_atoms_message(stats, resolved, popt), file=sys.stderr)
        return EXIT_NO_ATOMS

    say(f"got <{structure.natoms}> atoms.")
    for line in stats.lines(popt.legacy_elements):
        say(line)

    if args.info:
        _print_info(structure, stats, stream)
        return EXIT_OK

    ident = args.name or pov_identifier(args.output if args.output != "-" else "molecule")
    sopt = _scene_options_from(args, ident)

    say("Re-orienting and positioning structure.")
    prepare_structure(structure, sopt)

    bonds = None
    if sopt.ball_stick:
        bonds = find_bonds(
            structure, sopt.bond_threshold, sopt.bond_mode, sopt.bond_tolerance
        )
        say(f"Computing bonds: found {len(bonds)} bonds.")
        if not bonds:
            print(
                "pdb2pov: no bonds found; try a larger -d threshold"
                if sopt.bond_mode is BondMode.DISTANCE
                else "pdb2pov: no bonds found; try a larger --bond-tolerance",
                file=sys.stderr,
            )
            return EXIT_NO_BONDS

    if args.output == "-":
        write_scene(structure, sopt, sys.stdout, bonds)
        return EXIT_OK

    path = args.output
    if not path.lower().endswith((".pov", ".inc")):
        path += sopt.extension

    say(f"Writing output file <{path}>")
    try:
        write_scene(structure, sopt, path, bonds)
    except OSError as exc:
        print(f"pdb2pov: can't write '{path}': {exc.strerror or exc}", file=sys.stderr)
        return EXIT_CANT_WRITE

    return EXIT_OK


def _no_atoms_message(stats, resolved: str, popt: ParseOptions) -> str:
    if stats.skipped_chain:
        return f"pdb2pov: no atoms left after filtering to chain(s) '{popt.chains}'"
    if stats.skipped_model:
        return (
            f"pdb2pov: no atoms in model {popt.model}; the file has "
            f"{stats.models_seen or 'none'}"
        )
    if stats.skipped_altloc:
        return (
            "pdb2pov: every atom was skipped as an alternate conformation.\n"
            "         Column 17 should be blank or 'A'; check the record layout, "
            "or pass\n         --altloc first."
        )
    if stats.skipped_malformed:
        return (
            f"pdb2pov: {stats.skipped_malformed} coordinate record(s) in "
            f"<{resolved}> could not be parsed"
        )
    return f"pdb2pov: no atoms found in <{resolved}>"


def _print_info(structure, stats, stream) -> None:
    print(f"  format:   {structure.source_format}", file=stream)
    if structure.title:
        print(f"  title:    {structure.title}", file=stream)
    print(f"  atoms:    {structure.natoms}", file=stream)
    chains = structure.chains()
    shown = ", ".join(repr(c) for c in chains[:20]) + (", ..." if len(chains) > 20 else "")
    print(f"  chains:   {len(chains)} ({shown})", file=stream)
    if stats.models_seen:
        print(f"  models:   {stats.models_seen}", file=stream)
    counts = structure.element_counts()
    summary = ", ".join(f"{sym} {n}" for sym, n in list(counts.items())[:15])
    print(f"  elements: {summary}", file=stream)

    extents = structure.extents()
    print(
        f"  extent:   x {extents.xmin:.3f} .. {extents.xmax:.3f}, "
        f"y {extents.ymin:.3f} .. {extents.ymax:.3f}, "
        f"z {extents.zmin:.3f} .. {extents.zmax:.3f}  (padded)",
        file=stream,
    )


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
