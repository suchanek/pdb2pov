/*
   Object and texture declarations for some stock Atoms. For POV-Ray
   Version 2.0
   Author:
	Eric G. Suchanek, Ph.D.
*/

/* Must include:
#include "colors.inc"
#include "shapes.inc"
#include "textures.inc"
 and either atoms_cpk.inc or atoms_covalent.inc...
*/

#declare Atoms_Inc_Temp = version

#version 2.0

#declare Clear_Atom = texture {
   pigment {color red 1.0 green 1.0 blue 1.0 filter 0.95}
   finish {Shiny
    ambient 0.2
    diffuse 0.0
    reflection 0.05
    refraction 1.0
    ior 1.0
   }
}

#declare Atom_N = object {
	sphere {<0.0, 0.0, 0.0> N_RAD }
	texture {
         finish { Shiny }
         pigment {
          color Blue
          quick_color Blue
        }
   }
}

#declare Atom_C = object {
	sphere {<0.0, 0.0, 0.0> C_RAD }
	texture {
	 finish {Shiny }
	 pigment {
  	  color Green 
	  quick_color Green
        }
   }
}

#declare Atom_O = object {
	sphere {<0.0, 0.0, 0.0> O_RAD }
        texture {
	 finish {Shiny }
	 pigment {
  	  color Red 
	  quick_color Red
        }
    }
}


#declare Atom_S = object {
	sphere {<0.0, 0.0, 0.0> S_RAD }
        texture {
	 finish {Shiny }
	 pigment {
  	  color Yellow 
	  quick_color Yellow
        }
    }

}

#declare Atom_P = object {
	sphere {<0.0, 0.0, 0.0> P_RAD }
        texture {
	 finish {Shiny }
	 pigment {
  	  color Yellow 
	  quick_color Yellow
        }
    }
}

#declare Atom_Fe = object {
	sphere {<0.0, 0.0, 0.0> FE_RAD }
        texture {
	 finish {Shiny }
	 pigment {
  	  color DarkPurple
	  quick_color DarkPurple
        }
    }
}

#declare Atom_Ca = object {
	sphere {<0.0, 0.0, 0.0> CA_RAD }
        texture {
	 finish {Shiny }
	 pigment {
  	  color White
	  quick_color White
        }
    }
}

#declare Atom_H = object {
	sphere {<0.0, 0.0, 0.0> H_RAD }
        texture {
	 finish {Shiny }
	 pigment {
  	  color VLightGrey
	  quick_color VLightGrey
        }
   }
}

#declare White_Bond = texture 
{
  finish {Shiny }
  pigment {
	   color VLightGrey
	   quick_color VLightGrey
          }
}



#version Atoms_Inc_Temp
