/*
 * atoms_glass2.inc -- glass-textured atoms, so surface atoms can be seen
 * through.  Used by pdb2pov's -q hybrid mode, which merges a glass shell over
 * a ball-and-stick model.
 *
 * Author: Eric G. Suchanek, Ph.D.
 * Copyright (c) 1993-2026 Eric G. Suchanek, Ph.D.  All rights reserved.
 *
 * Include after atoms2.inc; see the note there about colours.
 *
 * Calcium and iron are deliberately not the same colour here as in
 * atoms2.inc -- calcium is white solid but dark purple in glass, iron dark
 * purple solid but gold in glass.  That inconsistency dates from the 1994
 * files and is preserved so existing -q renders do not change.
 */

#declare Atoms_Glass_Inc_Temp = version;

#version 3.7;

#declare my_glass_finish = finish {
    ambient 0.1
    diffuse 0.1
    specular 0.8
    roughness 0.0003
    reflection 0.1
    phong 1
    phong_size 400
}

/* Attached to each glass sphere below; carries what finish{} used to. */
#declare my_glass_interior = interior { ior 1.15 }


#declare Glass_H = texture {
    pigment { color VLightGrey filter 0.98 }
    finish { my_glass_finish }
}

#declare Glass_C = texture {
    pigment { color rgbf <0.0, 1.0, 0.0, 0.98> }
    finish { my_glass_finish }
}

#declare Glass_N = texture {
    pigment { color rgbf <0.0, 0.0, 1.0, 0.98> }
    finish { my_glass_finish }
}

#declare Glass_O = texture {
    pigment { color rgbf <1.0, 0.0, 0.0, 0.98> }
    finish { my_glass_finish }
}

#declare Glass_S = texture {
    pigment { color rgbf <1.0, 1.0, 0.0, 0.98> }
    finish { my_glass_finish }
}

#declare Glass_P = texture {
    pigment { color rgbf <1.0, 1.0, 0.0, 0.98> }
    finish { my_glass_finish }
}

#declare Glass_Ca = texture {
    pigment { color DarkPurple filter 0.98 }
    finish { my_glass_finish }
}

#declare Glass_Fe = texture {
    pigment { color Gold filter 0.98 }
    finish { my_glass_finish }
}

#declare Glass_Se = texture {
    pigment { color rgbf <1.000, 0.631, 0.000, 0.98> }
    finish { my_glass_finish }
}

#declare Glass_F = texture {
    pigment { color rgbf <0.565, 0.878, 0.314, 0.98> }
    finish { my_glass_finish }
}

#declare Glass_Cl = texture {
    pigment { color rgbf <0.122, 0.941, 0.122, 0.98> }
    finish { my_glass_finish }
}

#declare Glass_Br = texture {
    pigment { color rgbf <0.651, 0.161, 0.161, 0.98> }
    finish { my_glass_finish }
}

#declare Glass_I = texture {
    pigment { color rgbf <0.580, 0.000, 0.580, 0.98> }
    finish { my_glass_finish }
}

#declare Glass_Li = texture {
    pigment { color rgbf <0.800, 0.502, 1.000, 0.98> }
    finish { my_glass_finish }
}

#declare Glass_Na = texture {
    pigment { color rgbf <0.671, 0.361, 0.949, 0.98> }
    finish { my_glass_finish }
}

#declare Glass_K = texture {
    pigment { color rgbf <0.561, 0.251, 0.831, 0.98> }
    finish { my_glass_finish }
}

#declare Glass_Mg = texture {
    pigment { color rgbf <0.541, 1.000, 0.000, 0.98> }
    finish { my_glass_finish }
}

#declare Glass_Mn = texture {
    pigment { color rgbf <0.611, 0.478, 0.780, 0.98> }
    finish { my_glass_finish }
}

#declare Glass_Co = texture {
    pigment { color rgbf <0.941, 0.565, 0.627, 0.98> }
    finish { my_glass_finish }
}

#declare Glass_Ni = texture {
    pigment { color rgbf <0.314, 0.816, 0.314, 0.98> }
    finish { my_glass_finish }
}

#declare Glass_Cu = texture {
    pigment { color rgbf <0.784, 0.502, 0.200, 0.98> }
    finish { my_glass_finish }
}

#declare Glass_Zn = texture {
    pigment { color rgbf <0.490, 0.502, 0.690, 0.98> }
    finish { my_glass_finish }
}

#declare Glass_Mo = texture {
    pigment { color rgbf <0.329, 0.710, 0.710, 0.98> }
    finish { my_glass_finish }
}

#declare Glass_W = texture {
    pigment { color rgbf <0.129, 0.580, 0.839, 0.98> }
    finish { my_glass_finish }
}

#declare Glass_Ag = texture {
    pigment { color rgbf <0.753, 0.753, 0.753, 0.98> }
    finish { my_glass_finish }
}

#declare Glass_Cd = texture {
    pigment { color rgbf <1.000, 0.851, 0.561, 0.98> }
    finish { my_glass_finish }
}

#declare Glass_Pt = texture {
    pigment { color rgbf <0.816, 0.816, 0.878, 0.98> }
    finish { my_glass_finish }
}

#declare Glass_Au = texture {
    pigment { color rgbf <1.000, 0.820, 0.137, 0.98> }
    finish { my_glass_finish }
}

#declare Glass_Hg = texture {
    pigment { color rgbf <0.722, 0.722, 0.816, 0.98> }
    finish { my_glass_finish }
}

#declare Glass_B = texture {
    pigment { color rgbf <1.000, 0.710, 0.710, 0.98> }
    finish { my_glass_finish }
}

#declare Glass_Si = texture {
    pigment { color rgbf <0.941, 0.784, 0.627, 0.98> }
    finish { my_glass_finish }
}

#declare Glass_As = texture {
    pigment { color rgbf <0.741, 0.502, 0.890, 0.98> }
    finish { my_glass_finish }
}

#declare Glass_Xe = texture {
    pigment { color rgbf <0.259, 0.620, 0.690, 0.98> }
    finish { my_glass_finish }
}


#declare Glass_X = texture {
    pigment { color Gray50 filter 0.98 }
    finish { my_glass_finish }
}

#declare Clear_Atom = texture {
    pigment { color rgbf <1.0, 1.0, 1.0, 0.98> }
    finish { my_glass_finish }
}


#declare Atom_Glass_H = sphere {
    <0, 0, 0>, H_RAD
    texture { Glass_H }
    interior { my_glass_interior }
}

#declare Atom_Glass_C = sphere {
    <0, 0, 0>, C_RAD
    texture { Glass_C }
    interior { my_glass_interior }
}

#declare Atom_Glass_N = sphere {
    <0, 0, 0>, N_RAD
    texture { Glass_N }
    interior { my_glass_interior }
}

#declare Atom_Glass_O = sphere {
    <0, 0, 0>, O_RAD
    texture { Glass_O }
    interior { my_glass_interior }
}

#declare Atom_Glass_S = sphere {
    <0, 0, 0>, S_RAD
    texture { Glass_S }
    interior { my_glass_interior }
}

#declare Atom_Glass_P = sphere {
    <0, 0, 0>, P_RAD
    texture { Glass_P }
    interior { my_glass_interior }
}

#declare Atom_Glass_Ca = sphere {
    <0, 0, 0>, CA_RAD
    texture { Glass_Ca }
    interior { my_glass_interior }
}

#declare Atom_Glass_Fe = sphere {
    <0, 0, 0>, FE_RAD
    texture { Glass_Fe }
    interior { my_glass_interior }
}

#declare Atom_Glass_Se = sphere {
    <0, 0, 0>, SE_RAD
    texture { Glass_Se }
    interior { my_glass_interior }
}

#declare Atom_Glass_F = sphere {
    <0, 0, 0>, F_RAD
    texture { Glass_F }
    interior { my_glass_interior }
}

#declare Atom_Glass_Cl = sphere {
    <0, 0, 0>, CL_RAD
    texture { Glass_Cl }
    interior { my_glass_interior }
}

#declare Atom_Glass_Br = sphere {
    <0, 0, 0>, BR_RAD
    texture { Glass_Br }
    interior { my_glass_interior }
}

#declare Atom_Glass_I = sphere {
    <0, 0, 0>, I_RAD
    texture { Glass_I }
    interior { my_glass_interior }
}

#declare Atom_Glass_Li = sphere {
    <0, 0, 0>, LI_RAD
    texture { Glass_Li }
    interior { my_glass_interior }
}

#declare Atom_Glass_Na = sphere {
    <0, 0, 0>, NA_RAD
    texture { Glass_Na }
    interior { my_glass_interior }
}

#declare Atom_Glass_K = sphere {
    <0, 0, 0>, K_RAD
    texture { Glass_K }
    interior { my_glass_interior }
}

#declare Atom_Glass_Mg = sphere {
    <0, 0, 0>, MG_RAD
    texture { Glass_Mg }
    interior { my_glass_interior }
}

#declare Atom_Glass_Mn = sphere {
    <0, 0, 0>, MN_RAD
    texture { Glass_Mn }
    interior { my_glass_interior }
}

#declare Atom_Glass_Co = sphere {
    <0, 0, 0>, CO_RAD
    texture { Glass_Co }
    interior { my_glass_interior }
}

#declare Atom_Glass_Ni = sphere {
    <0, 0, 0>, NI_RAD
    texture { Glass_Ni }
    interior { my_glass_interior }
}

#declare Atom_Glass_Cu = sphere {
    <0, 0, 0>, CU_RAD
    texture { Glass_Cu }
    interior { my_glass_interior }
}

#declare Atom_Glass_Zn = sphere {
    <0, 0, 0>, ZN_RAD
    texture { Glass_Zn }
    interior { my_glass_interior }
}

#declare Atom_Glass_Mo = sphere {
    <0, 0, 0>, MO_RAD
    texture { Glass_Mo }
    interior { my_glass_interior }
}

#declare Atom_Glass_W = sphere {
    <0, 0, 0>, W_RAD
    texture { Glass_W }
    interior { my_glass_interior }
}

#declare Atom_Glass_Ag = sphere {
    <0, 0, 0>, AG_RAD
    texture { Glass_Ag }
    interior { my_glass_interior }
}

#declare Atom_Glass_Cd = sphere {
    <0, 0, 0>, CD_RAD
    texture { Glass_Cd }
    interior { my_glass_interior }
}

#declare Atom_Glass_Pt = sphere {
    <0, 0, 0>, PT_RAD
    texture { Glass_Pt }
    interior { my_glass_interior }
}

#declare Atom_Glass_Au = sphere {
    <0, 0, 0>, AU_RAD
    texture { Glass_Au }
    interior { my_glass_interior }
}

#declare Atom_Glass_Hg = sphere {
    <0, 0, 0>, HG_RAD
    texture { Glass_Hg }
    interior { my_glass_interior }
}

#declare Atom_Glass_B = sphere {
    <0, 0, 0>, B_RAD
    texture { Glass_B }
    interior { my_glass_interior }
}

#declare Atom_Glass_Si = sphere {
    <0, 0, 0>, SI_RAD
    texture { Glass_Si }
    interior { my_glass_interior }
}

#declare Atom_Glass_As = sphere {
    <0, 0, 0>, AS_RAD
    texture { Glass_As }
    interior { my_glass_interior }
}

#declare Atom_Glass_Xe = sphere {
    <0, 0, 0>, XE_RAD
    texture { Glass_Xe }
    interior { my_glass_interior }
}


/* Generic atom; see the Atom_X note in atoms2.inc. */
#declare Atom_Glass_X = sphere {
    <0, 0, 0>, X_RAD
    texture { Glass_X }
    interior { my_glass_interior }
}

#version Atoms_Glass_Inc_Temp;
