/*
 * atoms_glass2.inc -- glass-textured atoms, so surface atoms can be seen
 * through.  Used by pdb2pov's -q hybrid mode, which merges a glass shell over
 * a ball-and-stick model.
 *
 * Author: Eric G. Suchanek, Ph.D.
 * Copyright (c) 1993-2026 Eric G. Suchanek, Ph.D.  All rights reserved.
 *
 * Include, in this order:
 *
 *   #include "colors.inc"
 *   #include "shapes.inc"
 *   #include "textures.inc"
 *   #include "atoms_vdw.inc"      // or atoms_covalent.inc / atoms_cpk.inc
 *   #include "atoms2.inc"
 *   #include "atoms_glass2.inc"
 *
 * Updated from POV-Ray 2.0 to 3.7 syntax:
 *   - `sphere { <0,0,0> RADIUS }` gained the comma 3.x requires.
 *   - float declarations gained terminating semicolons.
 *   - refraction and ior are no longer accepted inside finish{}.  The old
 *     `refraction 0.9 ior 1.15` becomes an interior{} attached to each
 *     sphere; the transparency itself is carried by the pigment's filter,
 *     which is what actually let light through even in 2.x.
 *   - Atom_Glass_Ca carried both a texture{} and a bare pigment{}.  In 3.x
 *     that is an error -- the object already has a texture -- so the
 *     redundant pigment, which duplicated Glass_Ca's colour anyway, is gone.
 *   - Glass_P is declared rather than borrowing Glass_S.  Phosphorus and
 *     sulphur are both yellow in atoms2.inc, so the appearance is unchanged.
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

#declare Glass_N = texture {
    pigment { color rgbf <0.0, 0.0, 1.0, 0.98> }
    finish { my_glass_finish }
}

#declare Glass_C = texture {
    pigment { color rgbf <0.0, 1.0, 0.0, 0.98> }
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

#declare Glass_H = texture {
    pigment { color VLightGrey filter 0.98 }
    finish { my_glass_finish }
}

#declare Glass_Fe = texture {
    pigment { color Gold filter 0.98 }
    finish { my_glass_finish }
}

#declare Clear_Atom = texture {
    pigment { color rgbf <1.0, 1.0, 1.0, 0.98> }
    finish { my_glass_finish }
}

#declare Atom_Glass_N = sphere {
    <0, 0, 0>, N_RAD
    texture { Glass_N }
    interior { my_glass_interior }
}

#declare Atom_Glass_C = sphere {
    <0, 0, 0>, C_RAD
    texture { Glass_C }
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

#declare Atom_Glass_H = sphere {
    <0, 0, 0>, H_RAD
    texture { Glass_H }
    interior { my_glass_interior }
}

#declare Atom_Glass_Fe = sphere {
    <0, 0, 0>, FE_RAD
    texture { Glass_Fe }
    interior { my_glass_interior }
}

#declare Glass_X = texture {
    pigment { color Gray50 filter 0.98 }
    finish { my_glass_finish }
}

/* Generic atom; see the Atom_X note in atoms2.inc. */
#declare Atom_Glass_X = sphere {
    <0, 0, 0>, X_RAD
    texture { Glass_X }
    interior { my_glass_interior }
}

#version Atoms_Glass_Inc_Temp;
