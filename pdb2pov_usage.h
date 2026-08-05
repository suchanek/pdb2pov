/*
 * pdb2pov_usage.h -- the usage message.
 *
 * Copyright (c) 1993-2026 Eric G. Suchanek, Ph.D.  All rights reserved.
 * Subject to the GNU License.
 */

#ifndef PDB2POV_USAGE_H
#define PDB2POV_USAGE_H

/* Two %s substitutions, both argv[0]. */
static const char *PDB2POV_USAGE =
    "\n"
    "pdb2pov " PDB2POV_VERSION " -- Brookhaven PDB to POV-Ray "
    POV_VERSION " scene conversion\n"
    "\n"
    "USAGE: %s InputFile OutputFile [options]\n"
    "\n"
    "  Filenames are given without extensions; .pdb (or .atm) and .pov\n"
    "  (or .inc) are appended automatically.\n"
    "\n"
    "  Representation:\n"
    "    -v          van der Waals radii (default)\n"
    "    -c          covalent radii\n"
    "    -b          ball and stick\n"
    "    -q          ball and stick with glass atoms\n"
    "    -d x.x      bond cutoff threshold in angstroms (default 2.2)\n"
    "    -r x.x      scale factor applied to all atomic radii\n"
    "\n"
    "  Scene:\n"
    "    -o          object only -- write a .inc with no camera or lights\n"
    "    -p          plain white sky, no ground\n"
    "    -s          cloudy sky\n"
    "    -g          plain ground\n"
    "    -h          checkered ground\n"
    "    -a          area light instead of a point light\n"
    "\n"
    "  Orientation:\n"
    "    -x deg      absolute rotation about X\n"
    "    -y deg      absolute rotation about Y\n"
    "    -z deg      absolute rotation about Z\n"
    "\n"
    "  Input:\n"
    "    -t          input is in .atm format rather than PDB\n"
    "\n"
    "  Example: %s crambin crambin -s -h -b -d 1.5 -x 90\n"
    "    Converts crambin.pdb to crambin.pov with a checkered ground and\n"
    "    sky, rotates 90 degrees about X, and renders ball and stick with a\n"
    "    1.5 angstrom bond cutoff.\n"
    "\n"
    "  Copyright (c) 1993-2026 Eric G. Suchanek, Ph.D.  All rights "
    "reserved.\n"
    "\n";

#endif /* PDB2POV_USAGE_H */
