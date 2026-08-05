/*
 * atoms_cpk.inc -- CPK radii, in angstroms.
 *
 * Include before atoms2.inc, which uses these to size its spheres.
 * POV-Ray 3.5 and later require the terminating semicolon on a float
 * declaration.
 */

#declare N_RAD  = 1.54;
#declare C_RAD  = 1.7;
#declare O_RAD  = 1.4;
#declare P_RAD  = 1.9;
#declare S_RAD  = 1.8;
#declare CA_RAD = 1.274;
#declare H_RAD  = 1.2;
#declare FE_RAD = 2.0;

/* Generic radius, for elements with no dedicated entry above. */
#declare X_RAD  = 1.7;
