"""
Option objects for parsing and for scene generation.

Two dataclasses rather than the C's single ``Options`` struct: reading a
structure and drawing one are independent, and a caller using this as a PDB
reader should not have to invent a bond threshold.
"""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum

__all__ = [
    "RadiiSet",
    "Backdrop",
    "Ground",
    "AltLocPolicy",
    "BondMode",
    "InputFormat",
    "ParseOptions",
    "SceneOptions",
]


class RadiiSet(str, Enum):
    """Which set of atomic radii the scene includes."""

    VDW = "vdw"
    COVALENT = "covalent"
    CPK = "cpk"

    @property
    def include_file(self) -> str:
        return {
            RadiiSet.VDW: "atoms_vdw.inc",
            RadiiSet.COVALENT: "atoms_covalent.inc",
            RadiiSet.CPK: "atoms_cpk.inc",
        }[self]


class Backdrop(str, Enum):
    """What, if anything, is behind the molecule."""

    NONE = "none"
    SKY = "sky"
    PLAIN = "plain"


class Ground(str, Enum):
    NONE = "none"
    PLAIN = "plain"
    CHECKER = "checker"


class AltLocPolicy(str, Enum):
    """
    How to handle a residue modelled in more than one conformation.

    Keeping every conformer puts overlapping spheres at nearly identical
    positions in the scene, plus spurious bonds between the A and B copies,
    so something has to be chosen.

    ``A``
        Keep the blank and ``A`` indicators.  What the C does, and the
        default here so the two agree.
    ``FIRST``
        Keep whichever indicator appears first for each atom.  Equivalent to
        ``A`` on conventionally ordered files, and correct on the ones that
        label conformers ``B``/``C`` with no ``A``.
    ``OCCUPANCY``
        Keep the highest-occupancy conformer of each atom, which is the
        crystallographer's answer to the question.  Ties go to the first.
    ``ALL``
        Keep everything, as pdb2pov did before 2.1.
    """

    A = "a"
    FIRST = "first"
    OCCUPANCY = "occupancy"
    ALL = "all"


class BondMode(str, Enum):
    """
    How bonds are decided.

    ``DISTANCE``
        Every pair closer than a single cutoff, the historical behaviour.
    ``COVALENT``
        Pairs closer than the sum of their covalent radii plus a tolerance,
        so a C-H at 1.1 A and an S-S at 2.05 A are both found without a
        cutoff that also bonds every carbon to its neighbour's neighbour.
    """

    DISTANCE = "distance"
    COVALENT = "covalent"


class InputFormat(str, Enum):
    AUTO = "auto"
    PDB = "pdb"
    CIF = "cif"
    ATM = "atm"


@dataclass
class ParseOptions:
    """Everything the readers need."""

    fmt: InputFormat = InputFormat.AUTO

    #: Chain IDs to keep, as a string of single characters (``"AB"``), or
    #: ``None`` for all.  mmCIF chain IDs can be longer than one character,
    #: so a list is accepted too.
    chains: str | list[str] | None = None

    #: Model number to keep, or ``None`` for the first model in the file.
    #: :data:`ALL_MODELS` keeps every model.
    model: int | None = None

    altloc: AltLocPolicy = AltLocPolicy.A

    keep_hetatm: bool = True
    keep_water: bool = True

    #: Guess elements from atom names in the pre-2.1 way, dropping anything
    #: the guess cannot place.  For reproducing an existing render.
    legacy_elements: bool = False

    #: Raise instead of skipping when a coordinate record will not parse.
    strict: bool = False

    #: Stop after this many accepted atoms.  ``None`` for no limit -- the C
    #: pre-counted records and sized a fixed array, which is why it had one.
    max_atoms: int | None = None

    def accepts_chain(self, chain_id: str) -> bool:
        if self.chains is None:
            return True
        if isinstance(self.chains, str):
            return chain_id in self.chains
        return chain_id in self.chains


#: Sentinel for :attr:`ParseOptions.model` meaning "every model".
ALL_MODELS = -1


@dataclass
class SceneOptions:
    """Everything the scene writer needs."""

    radii: RadiiSet = RadiiSet.VDW
    backdrop: Backdrop = Backdrop.NONE
    ground: Ground = Ground.NONE

    ball_stick: bool = False
    glass_atoms: bool = False
    area_light: bool = False
    object_only: bool = False

    radii_scale: float = 1.0
    bond_threshold: float = 2.2
    bond_mode: BondMode = BondMode.DISTANCE
    bond_tolerance: float = 0.45

    xrot: float = 0.0
    yrot: float = 0.0
    zrot: float = 0.0

    #: POV-Ray identifier stem for the declared objects.  Derived from the
    #: output path when the CLI builds this.
    name: str = "molecule"

    #: Omit atoms with no dedicated texture instead of drawing them as the
    #: generic ``Atom_X``.  Set with ``--legacy-elements``: reproducing an
    #: old render means reproducing the omission too.  The drop is counted
    #: and reported either way, which is the part that was missing before.
    legacy_drop_unknown: bool = False

    #: Set false to leave the timestamp out of the header, which is what
    #: makes two runs byte-comparable.
    timestamp: bool = True

    @property
    def extension(self) -> str:
        return ".inc" if self.object_only else ".pov"
