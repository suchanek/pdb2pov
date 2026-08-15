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
 * The first eight elements keep their 1994 colours, which are not the CPK
 * convention -- carbon is green, phosphorus yellow, iron dark purple, calcium
 * white.  Elements added in 2.2 use Jmol/CPK colours.  Changing the original
 * eight would alter every existing render, so they stay as they are.
 *
 * Every element here must also appear in the ELEMENTS table in pdb2pov.c and
 * in each radius include; make check renders one atom of each to prove it.
 */

#declare Atoms_Inc_Temp = version;

#version 3.7;

/*
 * A transparent atom.  Attach Clear_Atom_Interior to any object textured with
 * Clear_Atom -- a texture cannot carry an interior, so the two are separate.
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


#declare Atom_H = sphere {
    <0, 0, 0>, H_RAD
    texture {
        finish { Shiny }
        pigment { color VLightGrey quick_color VLightGrey }
    }
}

#declare Atom_C = sphere {
    <0, 0, 0>, C_RAD
    texture {
        finish { Shiny }
        pigment { color Green quick_color Green }
    }
}

#declare Atom_N = sphere {
    <0, 0, 0>, N_RAD
    texture {
        finish { Shiny }
        pigment { color Blue quick_color Blue }
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

#declare Atom_Ca = sphere {
    <0, 0, 0>, CA_RAD
    texture {
        finish { Shiny }
        pigment { color White quick_color White }
    }
}

#declare Atom_Fe = sphere {
    <0, 0, 0>, FE_RAD
    texture {
        finish { Shiny }
        pigment { color DarkPurple quick_color DarkPurple }
    }
}

#declare Atom_Se = sphere {   // 2.2
    <0, 0, 0>, SE_RAD
    texture {
        finish { Shiny }
        pigment { color rgb <1.000, 0.631, 0.000> quick_color rgb <1.000, 0.631, 0.000> }
    }
}

#declare Atom_F = sphere {   // 2.2
    <0, 0, 0>, F_RAD
    texture {
        finish { Shiny }
        pigment { color rgb <0.565, 0.878, 0.314> quick_color rgb <0.565, 0.878, 0.314> }
    }
}

#declare Atom_Cl = sphere {   // 2.2
    <0, 0, 0>, CL_RAD
    texture {
        finish { Shiny }
        pigment { color rgb <0.122, 0.941, 0.122> quick_color rgb <0.122, 0.941, 0.122> }
    }
}

#declare Atom_Br = sphere {   // 2.2
    <0, 0, 0>, BR_RAD
    texture {
        finish { Shiny }
        pigment { color rgb <0.651, 0.161, 0.161> quick_color rgb <0.651, 0.161, 0.161> }
    }
}

#declare Atom_I = sphere {   // 2.2
    <0, 0, 0>, I_RAD
    texture {
        finish { Shiny }
        pigment { color rgb <0.580, 0.000, 0.580> quick_color rgb <0.580, 0.000, 0.580> }
    }
}

#declare Atom_Li = sphere {   // 2.2
    <0, 0, 0>, LI_RAD
    texture {
        finish { Shiny }
        pigment { color rgb <0.800, 0.502, 1.000> quick_color rgb <0.800, 0.502, 1.000> }
    }
}

#declare Atom_Na = sphere {   // 2.2
    <0, 0, 0>, NA_RAD
    texture {
        finish { Shiny }
        pigment { color rgb <0.671, 0.361, 0.949> quick_color rgb <0.671, 0.361, 0.949> }
    }
}

#declare Atom_K = sphere {   // 2.2
    <0, 0, 0>, K_RAD
    texture {
        finish { Shiny }
        pigment { color rgb <0.561, 0.251, 0.831> quick_color rgb <0.561, 0.251, 0.831> }
    }
}

#declare Atom_Mg = sphere {   // 2.2
    <0, 0, 0>, MG_RAD
    texture {
        finish { Shiny }
        pigment { color rgb <0.541, 1.000, 0.000> quick_color rgb <0.541, 1.000, 0.000> }
    }
}

#declare Atom_Mn = sphere {   // 2.2
    <0, 0, 0>, MN_RAD
    texture {
        finish { Shiny }
        pigment { color rgb <0.611, 0.478, 0.780> quick_color rgb <0.611, 0.478, 0.780> }
    }
}

#declare Atom_Co = sphere {   // 2.2
    <0, 0, 0>, CO_RAD
    texture {
        finish { Shiny }
        pigment { color rgb <0.941, 0.565, 0.627> quick_color rgb <0.941, 0.565, 0.627> }
    }
}

#declare Atom_Ni = sphere {   // 2.2
    <0, 0, 0>, NI_RAD
    texture {
        finish { Shiny }
        pigment { color rgb <0.314, 0.816, 0.314> quick_color rgb <0.314, 0.816, 0.314> }
    }
}

#declare Atom_Cu = sphere {   // 2.2
    <0, 0, 0>, CU_RAD
    texture {
        finish { Shiny }
        pigment { color rgb <0.784, 0.502, 0.200> quick_color rgb <0.784, 0.502, 0.200> }
    }
}

#declare Atom_Zn = sphere {   // 2.2
    <0, 0, 0>, ZN_RAD
    texture {
        finish { Shiny }
        pigment { color rgb <0.490, 0.502, 0.690> quick_color rgb <0.490, 0.502, 0.690> }
    }
}

#declare Atom_Mo = sphere {   // 2.2
    <0, 0, 0>, MO_RAD
    texture {
        finish { Shiny }
        pigment { color rgb <0.329, 0.710, 0.710> quick_color rgb <0.329, 0.710, 0.710> }
    }
}

#declare Atom_W = sphere {   // 2.2
    <0, 0, 0>, W_RAD
    texture {
        finish { Shiny }
        pigment { color rgb <0.129, 0.580, 0.839> quick_color rgb <0.129, 0.580, 0.839> }
    }
}

#declare Atom_Ag = sphere {   // 2.2
    <0, 0, 0>, AG_RAD
    texture {
        finish { Shiny }
        pigment { color rgb <0.753, 0.753, 0.753> quick_color rgb <0.753, 0.753, 0.753> }
    }
}

#declare Atom_Cd = sphere {   // 2.2
    <0, 0, 0>, CD_RAD
    texture {
        finish { Shiny }
        pigment { color rgb <1.000, 0.851, 0.561> quick_color rgb <1.000, 0.851, 0.561> }
    }
}

#declare Atom_Pt = sphere {   // 2.2
    <0, 0, 0>, PT_RAD
    texture {
        finish { Shiny }
        pigment { color rgb <0.816, 0.816, 0.878> quick_color rgb <0.816, 0.816, 0.878> }
    }
}

#declare Atom_Au = sphere {   // 2.2
    <0, 0, 0>, AU_RAD
    texture {
        finish { Shiny }
        pigment { color rgb <1.000, 0.820, 0.137> quick_color rgb <1.000, 0.820, 0.137> }
    }
}

#declare Atom_Hg = sphere {   // 2.2
    <0, 0, 0>, HG_RAD
    texture {
        finish { Shiny }
        pigment { color rgb <0.722, 0.722, 0.816> quick_color rgb <0.722, 0.722, 0.816> }
    }
}

#declare Atom_B = sphere {   // 2.2
    <0, 0, 0>, B_RAD
    texture {
        finish { Shiny }
        pigment { color rgb <1.000, 0.710, 0.710> quick_color rgb <1.000, 0.710, 0.710> }
    }
}

#declare Atom_Si = sphere {   // 2.2
    <0, 0, 0>, SI_RAD
    texture {
        finish { Shiny }
        pigment { color rgb <0.941, 0.784, 0.627> quick_color rgb <0.941, 0.784, 0.627> }
    }
}

#declare Atom_As = sphere {   // 2.2
    <0, 0, 0>, AS_RAD
    texture {
        finish { Shiny }
        pigment { color rgb <0.741, 0.502, 0.890> quick_color rgb <0.741, 0.502, 0.890> }
    }
}

#declare Atom_Xe = sphere {   // 2.2
    <0, 0, 0>, XE_RAD
    texture {
        finish { Shiny }
        pigment { color rgb <0.259, 0.620, 0.690> quick_color rgb <0.259, 0.620, 0.690> }
    }
}


/*
 * Generic atom, for elements with no dedicated texture above.  pdb2pov routes
 * anything unrecognised here rather than dropping it, so an exotic heteroatom
 * shows up as a neutral grey sphere instead of silently vanishing.
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
