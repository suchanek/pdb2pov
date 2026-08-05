/*
 * atoms2.inc -- object and texture declarations for stock atoms.
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
 *
 * Updated from POV-Ray 2.0 to 3.7 syntax:
 *   - `sphere { <0,0,0> RADIUS }` gained the comma 3.x requires between the
 *     centre and the radius.
 *   - float declarations gained terminating semicolons.
 *   - Clear_Atom's refraction and ior moved out of finish{}, which no longer
 *     accepts them, into a companion interior{}.
 *   - the spheres are declared directly rather than wrapped in object{}.
 */

#declare Atoms_Inc_Temp = version;

#version 3.7;

/*
 * A transparent atom.  Attach Clear_Atom_Interior to any object textured with
 * Clear_Atom -- a texture cannot carry an interior, so the two are separate:
 *
 *   sphere { <0,0,0>, 1 texture { Clear_Atom } interior { Clear_Atom_Interior } }
 */
#declare Clear_Atom = texture {
    pigment { color rgbf <1.0, 1.0, 1.0, 0.95> }
    finish {
        Shiny
        ambient 0.2
        diffuse 0.0
        reflection 0.05
    }
}

#declare Clear_Atom_Interior = interior { ior 1.0 }

#declare Atom_N = sphere {
    <0, 0, 0>, N_RAD
    texture {
        finish { Shiny }
        pigment { color Blue quick_color Blue }
    }
}

#declare Atom_C = sphere {
    <0, 0, 0>, C_RAD
    texture {
        finish { Shiny }
        pigment { color Green quick_color Green }
    }
}

#declare Atom_O = sphere {
    <0, 0, 0>, O_RAD
    texture {
        finish { Shiny }
        pigment { color Red quick_color Red }
    }
}

#declare Atom_S = sphere {
    <0, 0, 0>, S_RAD
    texture {
        finish { Shiny }
        pigment { color Yellow quick_color Yellow }
    }
}

#declare Atom_P = sphere {
    <0, 0, 0>, P_RAD
    texture {
        finish { Shiny }
        pigment { color Yellow quick_color Yellow }
    }
}

#declare Atom_Fe = sphere {
    <0, 0, 0>, FE_RAD
    texture {
        finish { Shiny }
        pigment { color DarkPurple quick_color DarkPurple }
    }
}

#declare Atom_Ca = sphere {
    <0, 0, 0>, CA_RAD
    texture {
        finish { Shiny }
        pigment { color White quick_color White }
    }
}

#declare Atom_H = sphere {
    <0, 0, 0>, H_RAD
    texture {
        finish { Shiny }
        pigment { color VLightGrey quick_color VLightGrey }
    }
}

/*
 * Generic atom, for elements with no dedicated texture above.  pdb2pov 2.1
 * and later route anything unrecognised here rather than dropping it, so an
 * ion or an exotic heteroatom shows up as a neutral grey sphere instead of
 * silently vanishing from the scene.
 */
#declare Atom_X = sphere {
    <0, 0, 0>, X_RAD
    texture {
        finish { Shiny }
        pigment { color Gray50 quick_color Gray50 }
    }
}

/* Applied to the cylinders pdb2pov emits between bonded atoms. */
#declare White_Bond = texture {
    finish { Shiny }
    pigment { color VLightGrey quick_color VLightGrey }
}

#version Atoms_Inc_Temp;
