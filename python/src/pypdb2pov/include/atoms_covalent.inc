/*
 * atoms_covalent.inc -- covalent radii
 *
 * Author: Eric G. Suchanek, Ph.D.
 * Copyright (c) 1993-2026 Eric G. Suchanek, Ph.D.  All rights reserved.
 *
 * Covalent radii, in angstroms.
 *
 * The eight elements pdb2pov has always drawn keep their historical values
 * (Pauling and the CAChe software, 5/4/94 -egs-).  Elements added in 2.2 use
 * Cordero et al. (2008); where a metal has spin-state-dependent radii the
 * low-spin value is used.
 *
 * Include before atoms2.inc, which uses these to size its spheres.
 */


#declare H_RAD  = 0.35;
#declare C_RAD  = 0.77;
#declare N_RAD  = 0.74;
#declare O_RAD  = 0.73;
#declare S_RAD  = 1.04;
#declare P_RAD  = 1.1;
#declare CA_RAD = 1.74;
#declare FE_RAD = 1.17;
#declare SE_RAD = 1.2;
#declare F_RAD  = 0.57;
#declare CL_RAD = 1.02;
#declare BR_RAD = 1.2;
#declare I_RAD  = 1.39;
#declare LI_RAD = 1.28;
#declare NA_RAD = 1.66;
#declare K_RAD  = 2.03;
#declare MG_RAD = 1.41;
#declare MN_RAD = 1.61;
#declare CO_RAD = 1.26;
#declare NI_RAD = 1.24;
#declare CU_RAD = 1.32;
#declare ZN_RAD = 1.22;
#declare MO_RAD = 1.54;
#declare W_RAD  = 1.62;
#declare AG_RAD = 1.45;
#declare CD_RAD = 1.44;
#declare PT_RAD = 1.36;
#declare AU_RAD = 1.36;
#declare HG_RAD = 1.32;
#declare B_RAD  = 0.84;
#declare SI_RAD = 1.11;
#declare AS_RAD = 1.19;
#declare XE_RAD = 1.4;

/* Generic radius, for elements with no dedicated entry above. */
#declare X_RAD  = 1;
