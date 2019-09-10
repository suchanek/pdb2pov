/* $Id: pdb2pov_usage.h,v 1.6 1994/07/24 21:18:10 suchanek Exp $ */

/* $Log: pdb2pov_usage.h,v $
 * Revision 1.6  1994/07/24  21:18:10  suchanek
 * stupid errors...
 *
 * Revision 1.5  1994/07/24  21:17:30  suchanek
 * *** empty log message ***
 *
 * Revision 1.4  1994/07/24  21:16:43  suchanek
 * some rewording
 *
 * Revision 1.3  1994/07/24  21:07:21  suchanek
 * added a shorter usage string for CPK
 *
 * Revision 1.2  1994/07/18  21:33:01  eric
 * added -p flag
 *
 * Revision 1.1  1993/11/25  19:59:43  eric
 * Initial revision
 * */

static char  *PDB2POV_USAGE = "\nProgram: pdb2pov.xxx 1.19 1994/7/19 22:00 \n\n\
\nUSAGE: %s InputFile OutputFile \n\
\t[-o object_only]\n\
\t[-t (atm file format)]\n\
\t[-p (no sky or ground)]\n\
\t[-s (writes cloudy sky)]\n\
\t[-g (writes plain ground)]\n\
\t[-h (writes checkered ground)]\n\
\t[-a (create area light)]\n\
\t[-v (do van der Waals radii)]\n\
\t[-c (do covalent radii)]\n\
\t[-b (do ball_and_stick)]\n\
\t[-d x.x (bond cutoff threshold)]\n\
\t[-q (ball and stick + glass atoms)]\n\
\t[-x X-Axis rotation]\n\
\t[-y Y-Axis rotation]\n\
\t[-z Z-Axis rotation]\n\n\n\
\tExample: %s crambin crambin -s -h -b -d 1.5 -x 90\n\
\tConverts crambin.pdb to crambin.pov, writes checkered ground, sky,\n\
\trotates the protein 90 degrees in X and renders ball and stick mode\n\
\twith a distance cutoff of 1.5 angstroms.\n\n\
\tCopyright (c) 1993-2014 Eric G. Suchanek, Ph.D.\n\tAll rights reserved.\n\n";

static char  *PDB2POV_HELLO = "\nProgram: %s \n\n\
USAGE: %s InputFile OutputFile \n\n\
[-o object_only]\n\
[-t (atm file format)]\n\
[-p (no sky or ground)]\n\
[-s (writes cloudy sky)]\n\
[-g (writes plain ground)]\n\
[-h (writes checkered ground)]\n\
[-a (create area light)]\n\
[-v (do van der Waals radii)]\n\
[-c (do covalent radii)]\n\
[-b (do ball_and_stick)]\n\
[-d x.x (bond cutoff threshold)]\n\
[-q (ball and stick + glass atoms)]\n\
[-x X-Axis rotation]\n\
[-y Y-Axis rotation]\n\
[-z Z-Axis rotation]\n\n\
Example: %s crambin crambin -s -h -b -d 1.5 -x 90\n\
Converts crambin.pdb to crambin.pov, writes checkered ground, sky,\n\
rotates the protein 90 degrees in X and renders ball and stick mode\n\
with a distance cutoff of 1.5 angstroms.\n\n\
Copyright (c) 1993-2014 Eric G. Suchanek, Ph.D.\n\
All rights reserved.\n";

static char  *PDB2POV_HELLO_NEW = "\nProgram: %s \n\n\
This program is capable of converting Brookhaven Protein \n\
Databank atomic structure files to POV-Ray v2.x scene files. \n\
The program supports a wide variety of rendering options.  See the \n\
documentation for more details. \n\
Copyright (c) 1993-2014 Eric G. Suchanek, Ph.D.\n\
All rights reserved.\n";



