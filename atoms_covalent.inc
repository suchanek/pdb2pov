/*
 * atoms_covalent.inc -- covalent radii, in angstroms.
 *
 * From Pauling and the CAChe software.  5/4/94 -egs-
 *
 * Include before atoms2.inc, which uses these to size its spheres.
 * POV-Ray 3.5 and later require the terminating semicolon on a float
 * declaration.
 */

#declare N_RAD  = 0.74;
#declare C_RAD  = 0.77;
#declare O_RAD  = 0.73;
#declare P_RAD  = 1.10;
#declare S_RAD  = 1.04;
#declare CA_RAD = 1.74;
#declare H_RAD  = 0.35;
#declare FE_RAD = 1.17;

/* Generic radius, for elements with no dedicated entry above. */
#declare X_RAD  = 1.0;
