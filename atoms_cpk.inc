/*
 * atoms_cpk.inc -- CPK radii
 *
 * Author: Eric G. Suchanek, Ph.D.
 * Copyright (c) 1993-2026 Eric G. Suchanek, Ph.D.  All rights reserved.
 *
 * CPK radii, in angstroms.  Identical to atoms_vdw.inc: the Corey-Pauling-
 * Koltun space-filling radii are the van der Waals radii.  No command-line
 * flag selects this file; see the note in pdb2pov.h.
 *
 * Include before atoms2.inc, which uses these to size its spheres.
 */


#declare H_RAD  = 1.2;
#declare C_RAD  = 1.7;
#declare N_RAD  = 1.54;
#declare O_RAD  = 1.4;
#declare S_RAD  = 1.8;
#declare P_RAD  = 1.9;
#declare CA_RAD = 1.274;
#declare FE_RAD = 2;
#declare SE_RAD = 1.9;
#declare F_RAD  = 1.47;
#declare CL_RAD = 1.75;
#declare BR_RAD = 1.85;
#declare I_RAD  = 1.98;
#declare LI_RAD = 1.82;
#declare NA_RAD = 2.27;
#declare K_RAD  = 2.75;
#declare MG_RAD = 1.73;
#declare MN_RAD = 2.05;
#declare CO_RAD = 2;
#declare NI_RAD = 1.63;
#declare CU_RAD = 1.4;
#declare ZN_RAD = 1.39;
#declare MO_RAD = 2.1;
#declare W_RAD  = 2.1;
#declare AG_RAD = 1.72;
#declare CD_RAD = 1.58;
#declare PT_RAD = 1.75;
#declare AU_RAD = 1.66;
#declare HG_RAD = 1.55;
#declare B_RAD  = 1.92;
#declare SI_RAD = 2.1;
#declare AS_RAD = 1.85;
#declare XE_RAD = 2.16;

/* Generic radius, for elements with no dedicated entry above. */
#declare X_RAD  = 1.7;
