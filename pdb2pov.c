/****************************************************************************/
/*           $Id: pdb2pov.c,v 1.19 1994/07/18 21:32:26 eric Exp $       */
/*  A program to convert brookhaven PDB format atomic files into the POV    */
/*  format files needed for the POV ray-tracing program.	            */
/*  In addition to the files for the source program, this code is dependent */
/*  on the color file atoms.inc to render CPK style spheres.		    */
/*  Author: Eric G. Suchanek, Ph.D.      				    */
/*  Copyright 1993 Eric G. Suchanek, Ph.D. All Rights Reserved              */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *  */

/* $Log: pdb2pov.c,v $
 * Revision 1.19  1994/07/18  21:32:26  eric
 * added the -p option for plain sky, no ground
 *
 * Revision 1.18  1994/05/04  22:01:46  eric
 * changed bond radius to 0.17 rather than 0.3
 *
 * Revision 1.17  1994/03/15  18:06:28  eric
 * updated version string
 *
 * Revision 1.16  1994/01/22  22:23:48  eric
 * changed to use make_rotmat rather than inc_rotmat. Now rotations
 * are absolute, not relative
 *
 * Revision 1.15  1993/11/26  09:12:24  eric
 * added -a option for area lights
 *
 * Revision 1.14  1993/11/25  19:58:23  eric
 * added additional error checking and return values
 *
 * Revision 1.12  1993/11/23  18:50:32  eric
 * changed the default to NOT render the checkered ground, and to
 * generate covalent radii
 *
 * Revision 1.11  1993/11/17  17:19:08  eric
 * added code to handle the ball&stick/glass atom hybrid
 *
 * Revision 1.8  1993/11/16  18:35:01  eric
 * many additions, ball and stick mode, conditional defines for UNIX/AMIGA, etc
 * */

#define FNSIZE 256
#define FMSIZE 256
#define FESIZE 256


#ifndef AMIGA
#include <strings.h>
#include <ctype.h>
#include <math.h>

#else

#include <exec/types.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#ifdef _M68881
#include <m68881.h>
#endif
#define  __USE_SYSBASE
#include <proto/exec.h>
#include <proto/all.h>

#endif 

#include <stdio.h>
#include <time.h>
#include <stdlib.h>

#include "pdb2light.h"
#include "pdb2pov_errors.h"

#ifdef AMIGA
#include "asyncio.h"
#include "asyncio_protos.h"
#include "pdb2pov_protos.h"
#include "util_protos.h"
#define DEBUG 1

long __stack = 8192;
#endif

static char rcsid[] = "$Id: pdb2pov.c,v 1.19 1994/07/18 21:32:26 eric Exp $";

static char vers[] = "\0$VER: pdb2pov 1.17 (3/15/94)";

#include "pdb2pov_usage.h"

/*
\t[-m fractional-margin] \n\
\t[-s image-size (0 or 1)]\n\
\t[-[xyz] number-of-frames]\n\
*/


#ifndef AMIGA
#define min(a,b) ((a) < (b) ? (a) : (b))
#define max(a,b) ((a) > (b) ? (a) : (b))

void compute_extents(int natoms, double **atoms, double *xmin, double *xmax, double *ymin, 
		double *ymax, double *zmin, double *zmax);

void make_rotmat(double rotmat[3][3], double gamma, double beta, double theta); 

void matmul(double new[3][3], double old[3][3]);

#endif

double xrot = 0.0;
double yrot = 0.0;
double zrot = 0.0;

#define BOND_RAD .17


void print_coords(int natoms, double **atoms)
{
  int i;

  for (i = 0; i < natoms; i++) 
    fprintf(stderr,"Atom %d: %6.2f %6.2f %6.2f\n",i,atoms[i][0],atoms[i][1],
	atoms[i][2]);
  fprintf(stderr,"\n");
}


/* function repositions y coordinates such that the origin (0,0) is in the */
/* lower left corner of the screen					   */
 
void flip_yaxis(int natoms, double **atoms, int maximum_y)
{
  register int index;

  for (index = 0; index < natoms; index++) 
    atoms[index][Y] = maximum_y -  atoms[index][Y];
}

/* correct for a left-handed display system */

void flip_zaxis(int natoms,double **atoms)
{
  register int index;

  for (index = 0; index < natoms; index++) 
    atoms[index][Z] = -atoms[index][Z];
}


int parse_inputs(argc,argv,in,out,inc,axis,frame_counter,increment,
	     radii_scale, obj_only, do_pdb, do_sky, do_ground, do_check,
	     do_vdw, do_cov, do_cpk, do_ball, threshold, do_hybrid, do_area, do_plain)
     int argc;
     char **argv;
     char *in;
     char *out;
     double *inc;
     int *axis;
     int *frame_counter;
     int *increment;
     double *radii_scale;
     int *obj_only;
     int *do_pdb;
     int *do_sky;
     int *do_ground;
     int *do_check;
     int *do_vdw;
     int *do_cov;
     int *do_cpk;
     int *do_ball;
     int *do_area;
     int *do_plain;
     double *threshold;
     int *do_hybrid;
{
  int i,frames,rot_increment;
  double fudge,rad_scale;
  double thresh = 0.0;

  if (argc < 3) {
    fprintf(stderr,PDB2POV_USAGE,argv[0],argv[0]);
    return(0);
  }

  *axis = X; frames = 1; fudge = 1.0; rad_scale = 1.0;
  *radii_scale = 1.0;

  for (i = 3; i < argc && (*argv[i] == '-'); ++i) {
	switch(*(argv[i]+1)) {
	case 'v':
	  *do_vdw = 1;
	  *do_cpk = 0;
	  *do_cov = 0;
	  break;
	case 'c':
	  *do_vdw = 0;
	  *do_cpk = 0;
	  *do_cov = 1;
	  break;
	case 'b':
	  *do_ball = 1;
	  break;
	case 'd':
	  sscanf(argv[++i], " %lf",&thresh);
	  *threshold = thresh;
	  break;
	case 'm':
	  fudge = atof(argv[++i]);
	  break;
	case 'x':
	case 'X':
	  sscanf(argv[++i]," %lf",&xrot);
	  break;
	case 'y':
	case 'Y':
	  sscanf(argv[++i]," %lf",&yrot);
	  break;
	case 'z':
	case 'Z':
	  sscanf(argv[++i]," %lf",&zrot);
	  break;
	case 'i':
	case 'I':
	  sscanf(argv[++i]," %d",&rot_increment);
	  break;
	case 'r':
	  rad_scale = atof(argv[++i]);
	  *radii_scale = rad_scale <= 0.0 ? 1.0 : rad_scale;
	  break; 
	case 'o':
	case 'O':
	  *obj_only = 1;
	  break;
	case 't':
	case 'T':
	  *do_pdb = 0;
	  break;
	case 'a':
	case 'A':
	  *do_area = 1;
	  break;
	case 'p':
	case 'P':
	  *do_check = 0;
	  *do_ground = 0;
	  *do_sky = 0;
	  *do_plain = 1;
	  break;
	case 's':
	case 'S':
	  *do_sky = 1;
	  break;
	  
	case 'g':
	case 'G':
	  *do_ground = 1;
	  *do_check = 0;
	  break;
	case 'h':
	case 'H':
	  *do_check = 1;
	  *do_ground = 0;
	  break;
	case 'q':
	case 'Q':
	  *do_ball = 1;
	  *do_hybrid = 1;
	  break;
	default:
	  fprintf(stderr,"Parse error: Argument %s not understood exiting....\n", argv[i]);
	  return(0);
	  break;
        }	/* switch */
      }		/* for    */ 
  
  if (*do_pdb == 1)
    {
      strcpy(in,argv[1]);
      strcat(in,".pdb");
    }
  else
    {
      strcpy(in,argv[1]);
      strcat(in,".atm");
    }

  strcpy(out,argv[2]);

  *frame_counter = frames;
  *increment = rot_increment;
  *inc = fudge <= 0.0 ? 1.0 : fudge;
  return(1);
}

/* function translates protein to LOCAL geometric center of coordinates */
 
void center_protein(int natoms, double **atoms, double center[])
{
  register int i;
  register double xcenter, ycenter, zcenter;

  xcenter = ycenter = zcenter = 0.0;

  for (i = 0; i < natoms; i++) {
    xcenter += atoms[i][0];
    ycenter += atoms[i][1];
    zcenter += atoms[i][2];
  }
  
  xcenter /= natoms;
  ycenter /= natoms;
  zcenter /= natoms;

  for (i = 0; i < natoms; i++) {
    atoms[i][0] -= xcenter; 
    atoms[i][1] -= ycenter;
    atoms[i][2] -= zcenter;
  }
  center[X] = xcenter; center[Y] = ycenter; center[Z] = zcenter;
}

void protein_cofmass(int natoms, double **atoms,double center[])
{

  register int i;
  register double xcenter, ycenter, zcenter;

  xcenter = ycenter = zcenter = 0.0;

  for (i = 0; i < natoms; i++) {
    xcenter += atoms[i][0];
    ycenter += atoms[i][1];
    zcenter += atoms[i][2];
  }
  
  xcenter /= natoms;
  ycenter /= natoms;
  zcenter /= natoms;

  center[X] = xcenter; center[Y] = ycenter; center[Z] = zcenter;
}

/* function translates protein (in COFmass coords) to the center of the
   screen (640,512) and then auto-scales to fill window */

 
void translate_protein(int natoms, double **atoms, int MAX_X, int MAX_Y)
{

  register int i;

  for (i = 0; i < natoms; i++) {
    atoms[i][0] += CENTER_X;
    atoms[i][1] += CENTER_Y;
  }
}


/* not really used anymore, but I hate to pitch algorithms :) */

double compute_max_scale(int natoms, double **atoms, int MAX_X, int MAX_Y)
{
  double xmax,xmin,ymax,ymin,zmax,zmin;
  double xextent,yextent,extent,max_scale;

  compute_extents(natoms,atoms,&xmin,&xmax,&ymin,&ymax,&zmin,&zmax);

  xextent = xmax - xmin;
  yextent = ymax - ymin;

  extent = xextent > yextent ? xextent : yextent;
  max_scale = xextent > yextent ? MAX_X * FACTOR / xextent : MAX_Y * FACTOR / yextent;
  max_scale = max_scale <= 0.0 ? 1.0 : max_scale;
  fprintf(stderr,"\nScale factor: %lf",max_scale);
  return(max_scale);
}

/* computes the bounding box for the atom list */
void compute_extents(int natoms, double **atoms, double *xmin, double *xmax, double *ymin, 
		double *ymax, double *zmin, double *zmax)
{
  double xm,xma,ym,yma,zm,zma;
  int i;

  xm = ym = zm = 9999.9;
  xma = yma = zma = -9999.9;

  for (i = 0; i < natoms; i++) 
    {
      xm = min(atoms[i][0],xm);
      xma = max(atoms[i][0],xma);
      ym = min(atoms[i][1],ym);
      yma = max(atoms[i][1],yma);
      zm = min(atoms[i][2],zm);
      zma = max(atoms[i][2],zma);
    }

  /* NOTE: MAX_RAD is the maximum atom radius in Angstroms ! */

  *xmin = xm - MAX_RAD; *xmax = xma + MAX_RAD;
  *ymin = ym - MAX_RAD; *ymax = yma + MAX_RAD;
  *zmin = zm - MAX_RAD; *zmax = zma + MAX_RAD;

#ifdef DBG
  fprintf(stderr,"\nXmin:%6.2f Xmax: %6.2f\nYmin: %6.2f Ymax: %6.2f\nZmin: %6.2f Zmax: %6.2f\n",*xmin,*xmax,*ymin,*ymax,*zmin,*zmax); 
#endif
}

/* function rotates an atom list about a specific axis and center point by 'amount'
   degrees. */

void rotate_protein(int natoms, double **atoms, double **patoms, 
		    double xrot, double yrot, double zrot)
{
   int i;
   double matrix[3][3] = {{1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, 
			        {0.0, 0.0, 1.0}};
   make_rotmat(matrix,xrot, yrot, zrot);

   for (i = 0; i < natoms; i++) 
      {
	patoms[i][0] = atoms[i][0];
	patoms[i][1] = atoms[i][1];
	patoms[i][2] = atoms[i][2];
      }

   for (i = 0; i < natoms; i++) 
     {
       double tmp[3];

       tmp[0] = patoms[i][0]; tmp[1] = patoms[i][1]; tmp[2] = patoms[i][2]; 

       matrix_times_vector(matrix, &tmp[0], &tmp[1], &tmp[2]);
       patoms[i][0] = tmp[0]; patoms[i][1] = tmp[1]; patoms[i][2] = tmp[2];
     }
 }

/* write the camera definitions out to the output file */
 
void write_camera(FILE *f, double xmin, double xmax, double ymin, 
		  double ymax, double zmin, double zmax)
{
  double xavg = (xmax - xmin) / 2.0;
  double yavg = (ymax - ymin) / 2.0;
  double zavg = (zmax - zmin) / 2.0;

  double theta = 22.0 / 57.298;
  double dist, distx, disty, distz;

  distx = -xavg / (tan(theta));
  disty = -yavg / (tan(theta));
  distz = -zavg / (tan(theta));

  dist = min(distx, disty);
  dist = min(dist, distz);

  fprintf(f,"camera {\n");
  fprintf(f,"   location < 0, 0, %.3lf > \n", dist);

  fprintf(f,"   direction < 0, 0, 1> \n");
  fprintf(f, "   up <0, 1, 0> \n");
  fprintf(f, "   right <4/3, 0, 0> \n");
  fprintf(f, "   look_at <0, 0, 0> \n");
  fprintf(f,"}\n");

}

/* write the lighting definitions out to the output file */

void write_light(FILE *f, double xmin, double xmax, double ymin, 
		 double ymax, double zmin, double zmax, int do_area)
{
  double xavg = (xmax - xmin) / 2.0;

  double theta = 20.0 / 57.298;
  double dist = 1.0;

  dist = -xavg / (tan(theta));


  fprintf(f,"object { \n");
  fprintf(f,"  light_source {\n");
  fprintf(f,"     <%.3lf, %.3lf, %.3lf> \n", xmax, ymax, dist);
  fprintf(f,"     color White \n");
  if (do_area)
    {
      fprintf(f,"     area_light <%.3lf, 0, 0>, <0, 0, %.3lf>, 9, 9 \n", xmax / 2.0, xmax / 2.0);
      fprintf(f,"     adaptive 1 \n");
      fprintf(f,"     jitter \n");
   }

  fprintf(f,"  }\n");
  fprintf(f,"}\n");

}

void write_sky(FILE *f)
{
  fprintf(f,"/* A nice gradient shaded blue sky with white coulds */ \n");
  fprintf(f,"sphere { <0, 0, 0>, 3000 \n");
  fprintf(f,"   pigment { \n");
  fprintf(f,"      gradient <0, 1, 0> \n");
  fprintf(f,"      colour_map { \n");
  fprintf(f,"         [0, 0.8  colour red .3   green 0.3 blue 1 \n");
  fprintf(f,"                  colour red 0.7 green 0.7 blue 1] \n");
  fprintf(f,"         [0.8, 1  colour red 0.7 green 0.7 blue 1 \n");
  fprintf(f,"                  colour red 0.90 green 0.9 blue 1] \n");
  fprintf(f,"      } \n");
  fprintf(f,"      scale <3000,  3000,  3000> \n");
  fprintf(f,"      quick_colour red 0.7  green 0.7 blue 1.0 \n");
  fprintf(f,"   } \n");
  fprintf(f,"   finish { \n");
  fprintf(f,"      ambient 0.7 \n");
  fprintf(f,"      diffuse 0   /* we don't want clouds casting shadows on the sky */ \n");
  fprintf(f,"   } \n");
  fprintf(f,"} \n");
  fprintf(f," \n");
  fprintf(f,"sphere { <0, 0, 0>, 2590 \n");
  fprintf(f,"   pigment { \n");
  fprintf(f,"      bozo \n");
  fprintf(f,"      turbulence 0.5 \n");
  fprintf(f,"      colour_map { \n");
  fprintf(f,"         [0,   0.6   colour Clear \n");
  fprintf(f,"                     colour Clear] \n");
  fprintf(f,"         [0.6, 0.8   colour Clear \n");
  fprintf(f,"                     colour White] \n");
  fprintf(f,"         [0.8, 1.001 colour White \n");
  fprintf(f,"                     colour red 0.8 green 0.8 blue 0.8] \n");
  fprintf(f,"      } \n");
  fprintf(f,"      quick_colour red 0.7 green 0.7 blue 1 \n");
  fprintf(f,"      scale <1000, 200, 1000> \n");
  fprintf(f,"   } \n");
  fprintf(f,"   finish {ambient 0.7 diffuse 0} \n");
  fprintf(f,"} \n");
}

void write_sky_plain(FILE *f)
{
  fprintf(f,"/* A plain white sky */ \n");
  fprintf(f,"sphere { <0, 0, 0>, 3000 \n");
  fprintf(f,"         pigment { colour White } \n");
  fprintf(f,"         finish { \n");
  fprintf(f,"           ambient 0.7 \n");
  fprintf(f,"           diffuse 0   \n");
  fprintf(f,"          } \n");
  fprintf(f,"       } \n");
}


/* write the ground definitions out to the output file */

void write_ground(FILE *f, double xmin, double xmax, double 
		  ymin, double ymax, double zmin, double zmax)
{
  double yavg = (ymax - ymin) / 2.0;

  fprintf(f,"object { \n");
  fprintf(f,"  plane { y, %3lf }\n",
	  -1.0 * (fabs(ymin) + fabs(yavg) / 2.0));
  fprintf(f," pigment {color RichBlue quick_color RichBlue} \n");
  fprintf(f," normal {bumps 0.2} \n");
  fprintf(f,"}\n");
}

void write_check(FILE *f, double xmin, double xmax, double 
		  ymin, double ymax, double zmin, double zmax)
{
  double xavg = (xmax - xmin) / 2.0;
  double yavg = (ymax - ymin) / 2.0;

  fprintf(f,"/* Create the beloved famous raytrace checkered floor */ \n");
  fprintf(f,"plane { y, %3lf \n", 	  -1.0 * (fabs(ymin) + fabs(yavg) / 2.0));
  fprintf(f,"      pigment { \n");
  fprintf(f,"      checker colour Black colour White \n");
  fprintf(f,"      scale %3lf \n", ((xavg + yavg) / 3.0));
  fprintf(f,"     } \n");
  fprintf(f,"     finish { \n");
  fprintf(f,"       ambient 0.2 \n");
  fprintf(f,"       diffuse 0.8 \n");
  fprintf(f,"     } \n");
  fprintf(f,"} \n");

}

/* compute the bounding sphere for the atom list. adds the constant MAX_RAD to
   account for the radii of the atoms (overestimate) */

double compute_sphere(int natoms, double **atoms)

{
  double rad = 0.0;
  double rad_max = -999.9;
  double x, y, z;
  int i;

  for (i = 0; i < natoms; i++)
    {
      x = atoms[i][0] ; y = atoms[i][1]; z = atoms[i][2];
      rad = sqrt((x * x + y * y + z * z)) + MAX_RAD;
      if (rad > rad_max) rad_max = rad;
    }


  return (rad_max);

}

int write_output(char *out, int natoms, double **atoms, int *atom_types, 
	     int **bonds, int nbonds, double rad_factor, int obj_only, 
	     int do_sky, int do_ground, int do_check,
	     int do_cpk, int do_cov, int do_vdw, int do_ball, int do_hybrid, int do_area, int do_plain)
{

  char drive[FNSIZE], path[FMSIZE], node[FNSIZE], ext[FESIZE];

  FILE *f;
  FILE *fopen();
  int i;
  char tmp[128];

  char numb[8];
  double xmin, ymin, xmax, ymax, zmin, zmax;
  double sphere_rad = 0.0;
  int ok = 1;

  long t;
  struct tm *gm;

  char date[10];
  char ttime[9];

  char time_string[20];

  time(&t);
  gm = localtime(&t);
  sprintf(date,"%02d/%02d/%02d \0", gm->tm_mon+1, gm->tm_mday, gm->tm_year);
  sprintf(ttime,"%02d:%02d:%02d\0", gm->tm_hour, gm->tm_min, gm->tm_sec);

  strcpy(time_string,date);
  strcat(time_string,ttime);

  if (obj_only)
    {
      sprintf(numb,".inc");
      strcpy(tmp,out);
      strcat(tmp,numb);
    }
  else
    {
      sprintf(numb,".pov");
      strcpy(tmp,out);
      strcat(tmp,numb);
    }

#ifdef AMIGA
  strsfn(out, drive, path, node, ext);
#endif

  compute_extents(natoms,atoms,&xmin,&xmax,&ymin,&ymax,&zmin,&zmax);  
  sphere_rad = compute_sphere(natoms,atoms);
  sphere_rad += .02 * sphere_rad; /* fudge factor for bounding sphere */

  if ((f = fopen(tmp,"w")) == NULL) {
	fprintf(stderr,"\nCan't open output file (%s) for writing...\n",out);
	return(ERR_CANT_WRITE_OUTPUT);
  }
  
  printf("Writing output file <%s> \n", tmp);

#ifdef AMIGA
  fprintf(f,"//\n// PDB2POV (Amiga) atom input prepared by pdb2pov %s \n",time_string);
#else
  fprintf(f,"//\n// PDB2POV (UNIX) atom input prepared by pdb2pov %s \n",time_string);
#endif

  fprintf(f,"// Author: Eric G. Suchanek, Ph.D. \n");
  fprintf(f,"//\tAtoms: %4d \n",natoms);
  fprintf(f,"//\tExtent:\tXmin: %.3f Xmax: %.3f, \n//\t\tYmin: %.3f, Ymax: %.3f \n", 
	  xmin, xmax, ymin, ymax);
  fprintf(f,"//\t\tZmin: %.3f Zmax: %.3f \n",zmin,zmax);
  fprintf(f,"//\tEnclosing Sphere: %.3f\n\n", sphere_rad);

  if (! obj_only)
    {
      fprintf(f,"#include \"colors.inc\"\n");
      fprintf(f,"#include \"shapes.inc\"\n");
      fprintf(f,"#include \"textures.inc\"\n");
      if (do_cpk)
	fprintf(f,"#include \"atoms_cpk.inc\"\n"); /* radii */
      else
	{
	  if (do_vdw)
	    fprintf(f,"#include \"atoms_vdw.inc\"\n"); /* radii */
	  else
	    {
	      if (do_cov)
		fprintf(f,"#include \"atoms_covalent.inc\"\n"); /* radii */
	    }
	}


      fprintf(f,"#include \"atoms2.inc\"\n"); /* atom objects */

      if (do_hybrid)
	fprintf(f,"#include \"atoms_glass2.inc\"\n"); /* atom objects */

      write_camera(f, xmin, xmax, ymin, ymax, zmin, zmax);
      write_light(f, xmin, xmax, ymin, ymax, zmin, zmax, do_area);
      if (do_sky)
	write_sky(f);
      if (do_plain)
	write_sky_plain(f);
      if (do_ground)
	write_ground(f, xmin, xmax, ymin, ymax, zmin, zmax);
      else
	if (do_check)
	  write_check(f, xmin, xmax, ymin, ymax, zmin, zmax);

    }

  if (do_ball)
    {
      fprintf(f, "#declare BOND_RADIUS = %.2lf \n", BOND_RAD);
    }

  if (do_hybrid)
    {
      fprintf(f, "#declare ATM_SCL_B = %.2lf \n", rad_factor / .3);
    }

  fprintf(f, "#declare ATM_SCL = %.2lf \n \n", rad_factor);



  /* do the glass atoms */
  if (do_hybrid)
    {
#ifdef AMIGA
      fprintf(f,"#declare %s_obj_glass = merge {\n",node); /* was out */
#else
      fprintf(f,"#declare %s_obj_glass = merge {\n",out); /* was out */
#endif      
      for (i = 0; i < natoms; i++)
	{
	  switch(atom_types[i])
	    {
	    case H_TYPE:
	      fprintf(f,"object { Atom_Glass_H ");
	      ok = 1;
	      break;
	    case N_TYPE:
	      fprintf(f,"object { Atom_Glass_N ");
	      ok = 1;
	      break;
	    case C_TYPE:
	      fprintf(f,"object { Atom_Glass_C ");
	      ok = 1;
	      break;
	    case P_TYPE:
	      fprintf(f,"object { Atom_Glass_P ");
	      ok = 1;
	      break;
	    case O_TYPE:
	      fprintf(f,"object { Atom_Glass_O ");
	      ok = 1;
	      break;
	    case S_TYPE:
	      fprintf(f,"object { Atom_Glass_S ");
	      ok = 1;
	      break;
	    case CA_TYPE:
	      fprintf(f,"object { Atom_Glass_Ca ");
	      ok = 1;
	      break;
	    case FE_TYPE:
	      fprintf(f,"object { Atom_Glass_Fe ");
	      ok = 1;
	      break;
	    default:
	      ok = 0;
	      break;
	    }  /* case */
	  
	  if (ok)
	    {
	      fprintf(f,"scale <ATM_SCL_B, ATM_SCL_B, ATM_SCL_B> ");
	      fprintf(f,"translate <%.3f, %.3f, %.3f> }\n",
		      atoms[i][0],atoms[i][1],atoms[i][2]);
	    }
	  
	} /* for */
      fprintf(f,"}\n");   /* merge */
      
    }

#ifdef AMIGA
  fprintf(f,"#declare %s_obj = union {\n",node);
#else
  fprintf(f,"#declare %s_obj = union {\n",out);
#endif

  /* do the atoms */
  for (i = 0; i < natoms; i++)
    {
      switch(atom_types[i])
	{
	  case H_TYPE:
	   fprintf(f,"object { Atom_H ");
	   ok = 1;
	   break;
	  case N_TYPE:
	   fprintf(f,"object { Atom_N ");
	   ok = 1;
	   break;
	  case C_TYPE:
	   fprintf(f,"object { Atom_C ");
	   ok = 1;
	   break;
	  case P_TYPE:
	   fprintf(f,"object { Atom_P ");
	   ok = 1;
	   break;
	  case O_TYPE:
	   fprintf(f,"object { Atom_O ");
	   ok = 1;
	   break;
	  case S_TYPE:
	   fprintf(f,"object { Atom_S ");
	   ok = 1;
	   break;
	  case CA_TYPE:
	   fprintf(f,"object { Atom_Ca ");
	   ok = 1;
	   break;
	 case FE_TYPE:
	   fprintf(f,"object { Atom_Fe ");
	   ok = 1;
	   break;
 	  default:
	   ok = 0;
	   break;
	 }  /* case */

      if (ok)
	{
	  fprintf(f,"scale <ATM_SCL, ATM_SCL, ATM_SCL> ");
	  fprintf(f,"translate <%.3f, %.3f, %.3f> }\n",
		  atoms[i][0],atoms[i][1],atoms[i][2]);
	}

    }


  /* do the bonds */
  if (do_ball)
    for (i = 0; i < nbonds; i++)
      {
	fprintf(f, "cylinder { <%.2lf, %.2lf, %.2lf>, <%.2lf, %.2lf, %.2lf>, BOND_RADIUS \
texture {White_Bond} } \n", atoms[bonds[i][0]][0], 
	                    atoms[bonds[i][0]][1], 
	                    atoms[bonds[i][0]][2],
   	                    atoms[bonds[i][1]][0], 
   	                    atoms[bonds[i][1]][1], 
   	                    atoms[bonds[i][1]][2]);
/*	                    BOND_RAD); */

      }

  if (do_hybrid)
    {
#ifdef AMIGA
      fprintf(f, "object { %s_obj_glass} \n", node);
#else
      fprintf(f, "object { %s_obj_glass} \n", out);
#endif
    }
  fprintf(f,"}\n");   /* composite */

  /* now write a composite with bounding sphere */
#ifdef AMIGA
  fprintf(f,"\n#declare %s = object {\n\tobject { %s_obj } ", node, node);
#else
  fprintf(f,"\n#declare %s = object {\n\tobject { %s_obj } ", out, out);
#endif
  fprintf(f,"\n\tbounded_by {sphere {<0 0 0> %.3f} }\n     } \n",sphere_rad); 

  if (! obj_only)
#ifdef AMIGA
    fprintf(f,"\nobject { %s } \n", node);
#else
    fprintf(f,"\nobject { %s } \n", out);
#endif
}


/*********************** New Matrix functions *************************/

/* increment the rotation matrix by x, y, z */

void make_rotmat(double rotmat[3][3], double gamma, double beta, double theta) /* x-axis y-axis z-axis */
{

  double g = gamma / 57.29;
  double b = beta / 57.29;
  double t = theta / 57.29;

  register double sing = sin(g);
  register double cosg = cos(g);
  register double sinb = sin(b);
  register double cosb = cos(b);
  register double sint = sin(t);
  register double cost = cos(t);

  rotmat[0][0] = cost * cosb; 
  rotmat[0][1] = sint * cosg + cost * sinb * sing;
  rotmat[0][2] = sint * sing - cost * sinb * cosg;
  rotmat[1][0] = -sint * cosb;
  rotmat[1][1] = cost * cosg - sint * sinb * sing;
  rotmat[1][2] = cost * sing + sint * sinb * cosg;
  rotmat[2][0] = sinb;
  rotmat[2][1] = -cosb * sing;
  rotmat[2][2] = cosb * cosg;
  
}

void inc_rotmat(double rotmat[3][3], double gamma, double beta, 
		double theta) /* x-axis y-axis z-axis */
{

  double g = gamma / 57.29;
  double b = beta / 57.29;
  double t = theta / 57.29;

  double newmat[3][3];

  register double sing = sin(g);
  register double cosg = cos(g);
  register double sinb = sin(b);
  register double cosb = cos(b);
  register double sint = sin(t);
  register double cost = cos(t);

  newmat[0][0] = cost * cosb; 
  newmat[0][1] = sint * cosg + cost * sinb * sing;
  newmat[0][2] = sint * sing - cost * sinb * cosg;
  newmat[1][0] = -sint * cosb;
  newmat[1][1] = cost * cosg - sint * sinb * sing;
  newmat[1][2] = cost * sing + sint * sinb * cosg;
  newmat[2][0] = sinb;
  newmat[2][1] = -cosb * sing;
  newmat[2][2] = cosb * cosg;
  
  matmul(newmat, rotmat);
}

/* multiply the old and new matrices; replace old with product */

void matmul(double new[3][3], double old[3][3])
{
  double tmp[3][3];

    tmp[0][0]= new[0][0]*old[0][0]+new[0][1]*old[1][0]+new[0][2]*old[2][0]; 
    tmp[0][1]= new[0][0]*old[0][1]+new[0][1]*old[1][1]+new[0][2]*old[2][1]; 
    tmp[0][2]= new[0][0]*old[0][2]+new[0][1]*old[1][2]+new[0][2]*old[2][2]; 

    tmp[1][0]= new[1][0]*old[0][0]+new[1][1]*old[1][0]+new[1][2]*old[2][0]; 
    tmp[1][1]= new[1][0]*old[0][1]+new[1][1]*old[1][1]+new[1][2]*old[2][1]; 
    tmp[1][2]= new[1][0]*old[0][2]+new[1][1]*old[1][2]+new[1][2]*old[2][2]; 

    tmp[2][0]= new[2][0]*old[0][0]+new[2][1]*old[1][0]+new[2][2]*old[2][0]; 
    tmp[2][1]= new[2][0]*old[0][1]+new[2][1]*old[1][1]+new[2][2]*old[2][1]; 
    tmp[2][2]= new[2][0]*old[0][2]+new[2][1]*old[1][2]+new[2][2]*old[2][2]; 

    old[0][0]= tmp[0][0];
    old[0][1]= tmp[0][1];
    old[0][2]= tmp[0][2];

    old[1][0]= tmp[1][0];
    old[1][1]= tmp[1][1];
    old[1][2]= tmp[1][2];

    old[2][0]= tmp[2][0];
    old[2][1]= tmp[2][1];
    old[2][2]= tmp[2][2];
}

int matrix_times_vector(double matrix[3][3], double *x, double *y, double *z)
{
  register double xx, yy, zz;

  xx = matrix[0][0] * *x + matrix[0][1] * *y + matrix[0][2] * *z;
  yy = matrix[1][0] * *x + matrix[1][1] * *y + matrix[1][2] * *z;
  zz = matrix[2][0] * *x + matrix[2][1] * *y + matrix[2][2] * *z;
  *x = xx; *y = yy; *z = zz;
  return(1);
}

#ifdef AMIGA
#define CRET 0x0D
#define LF 0x0A

int getline(char *linebuf, int numb, struct AsyncFile *file)
{
  int i = 0;
  LONG read = 0;

  while (((read = ReadCharAsync(file)) >= 0) && (i < numb))
    {
      if (read == LF || read == CRET || read == 0)
	{
	  if (read == CRET)
	    read = ReadCharAsync(file); /* eat the linefeed */
	  if (read < 0)
	    return(0);

	  linebuf[i] = '\0';
	  return(i);
	}
      else
	{
	  if (read < 0)
	    {
	      return(NULL);
	    }
	  else
	    linebuf[i++] = (char)read;
	}
    }
}
#endif

/*  PDB file I/O */

int getpdb(char *fname, int maxnatoms, double **apos, 
	   char **atype, double *charge, char **amino)
{
  int natoms;
#ifdef AMIGA
  struct AsyncFile *file_in;
  LONG IO_buffersize = 32768;
  
  if (!(file_in = OpenAsync(fname, MODE_READ, IO_buffersize)))
    {
      fprintf(stderr,"can't open file %s\n",fname);
      return(-1);
    }
#else
  FILE *file_in;
  FILE *fopen();

  file_in = fopen(fname, "r");
  if(file_in==NULL) {
    fprintf(stderr,"can't open file %s\n",fname);
    return -1;
  }

#endif

  natoms = fgetpdb(file_in, maxnatoms, apos, atype, charge,amino);
  if (natoms > maxnatoms)
    {
      fprintf(stderr,"Too many atoms! Truncated to <%d> \n", maxnatoms);
    }

  if (file_in)
#ifdef AMIGA
    CloseAsync(file_in);
#else
  fclose(file_in);
#endif

  return natoms;
}

#ifdef AMIGA
int fgetpdb(struct AsyncFile *in, int maxnatoms, 
	    double **apos, char **atype, double *charge, char **amino)
{
#else
int fgetpdb(FILE *f, int maxnatoms, 
	    double **apos, char **atype, double *charge, char **amino)
{
#endif
/*
 * Read a Brookhaven file putting the coords into apos and the char
 *  string names into atype. Return number of atoms.  M Pique
 * Stop reading on End-of-file, or 'END' line, or line with just a ".".
 * Here is some sample input: 
ATOM     19  HN  THR     2      15.386  10.901   4.600
ATOM     20  HOG THR     2      14.161   9.481   4.420
ATOM     21  N   CYS     3      13.507  11.238   8.398
ATOM     22  CA  CYS     3      13.659  10.715   9.763
ATOM     23  C   CYS     3      12.265  10.443  10.307
HETATM   99  FE  XXX     4      33.265  10.443  10.307
ATOM         OT  OME    52     -16.009  30.871 -13.037                -0.543
HETATM       CU+2CU    152     -20.577   2.601 -14.587                 2.000
HETATM       ZN+2ZN    152     -17.911   1.372 -19.974                 2.000
END
                                                         These charge values
							 are non-standard
							 and are not genuine
							 PDB. They are optional.
 */
  char linebuf[200];
  register int n; /* number of atoms */
  double tcharge;

  n=0;
#ifdef AMIGA
  while(getline(linebuf, sizeof(linebuf), in)!=NULL && n < maxnatoms &&
#else
  while(fgets(linebuf, sizeof linebuf, f)!=NULL && n < maxnatoms &&
#endif
	 0!=strcmp(".", linebuf) &&
	 0!=strncmp("END", linebuf,3) &&
	 0!=strncmp("end", linebuf,3)
	 ) {
    printf("%s\n", linebuf);
		if( (0==strncmp("ATOM",linebuf,4)||
			0==strncmp("atom",linebuf,4) ||
			0==strncmp("HETATM",linebuf,6) ||
			0==strncmp("hetatm",linebuf,6) )
		&& 1==sscanf(&linebuf[12]," %3s", atype[n])
		&& 1==sscanf(&linebuf[17],"%3s", amino[n])
		&& 3==sscanf(&linebuf[30],"%lf %lf %lf",
		&apos[n][0],  &apos[n][1], &apos[n][2]) ) {
#ifdef DEBUG
			printf("atom %d type %s %s at %lf %lf %lf\n",
			  n, atype[n], amino[n], apos[n][0],  apos[n][1], apos[n][2]);
#endif
			if(1==sscanf(&linebuf[70],"%lf", &tcharge))
			 charge[n] = tcharge;
			else charge[n] = 0.00;
			if ((0 == strncmp("ATOM", linebuf,4)) && ((0 == strncmp(atype[n],"ca",2))
								  || (0 == strncmp(atype[n],"CA",2))))
			  atype[n][1] = ' ';
			n++;
			}
		}
	return n;
	}

#ifdef AMIGA
int fgetpdbxyz(struct AsyncFile *in, int maxnatoms, double **apos)
#else
int fgetpdbxyz(FILE *in, int maxnatoms, double **apos)
#endif
{
/*
 * Read a Brookhaven file putting the coords into apos only
 * Return number of atoms read before end-of-file or line with END or .
 */
char linebuf[200];
register int n; /* number of atoms */

	n=0;

#ifdef AMIGA
  while(getline(linebuf, sizeof(linebuf), in)!=NULL && n < maxnatoms &&
#else
  while(fgets(linebuf, sizeof linebuf, in)!=NULL && n < maxnatoms &&
#endif
	 0!=strcmp(".", linebuf) &&
	 0!=strncmp("END", linebuf,3) &&
	 0!=strncmp("end", linebuf,3)
	) {
		if( (0==strncmp("ATOM",linebuf,4)||
			0==strncmp("atom",linebuf,4) ||
			0==strncmp("HETATM",linebuf,6) ||
			0==strncmp("hetatm",linebuf,6) )
		&& 3==sscanf(&linebuf[30],"%lf %lf %lf",
		&apos[n][0],  &apos[n][1], &apos[n][2]) ) {
#ifdef DEBUG
			printf("atom %d at %lf %lf %lf\n",
			  n, apos[n][0],  apos[n][1], apos[n][2]);
#endif
			n++;
			}
		}
	return n;
	}


getpdb_atnum(char *fname)
{
  int natoms;
#ifdef AMIGA
  struct AsyncFile *file_in;
  LONG IO_buffersize = 32768;
  
  if (!(file_in = OpenAsync(fname, MODE_READ, IO_buffersize)))
    {
      fprintf(stderr,"can't open file %s\n",fname);
      return(ERR_CANT_READ_INPUT);
    }
#else

  FILE *file_in;

  file_in = fopen(fname, "r");
  if(file_in==NULL) {
    fprintf(stderr,"can't open file %s\n",fname);
    return ERR_CANT_READ_INPUT;
  }

#endif

  natoms = fgetpdb_atnum(file_in);

  if (file_in)
#ifdef AMIGA
    CloseAsync(file_in);
#else
    fclose(file_in);
#endif
  return natoms;
}

#ifdef AMIGA
int fgetpdb_atnum(struct AsyncFile *in)
#else
int fgetpdb_atnum(FILE *in)
#endif
{
  char linebuf[200];
  register int n; /* number of atoms */

  n=0;

#ifdef AMIGA
  while(getline(linebuf, sizeof(linebuf), in)!=NULL && 
#else
  while(fgets(linebuf, sizeof linebuf, in)!=NULL && 
#endif
	0!=strcmp(".", linebuf) && 0!=strncmp("END", linebuf,3) &&
	0!=strncmp("end", linebuf,3) )
    {
      if( (0==strncmp("ATOM",linebuf,4)||
	   0==strncmp("atom",linebuf,4) ||
	   0==strncmp("HETATM",linebuf,6) ||
	   0==strncmp("hetatm",linebuf,6) ) )
	{
	  n++;
	}
    }
  return n;
}


/* .atm File I/O */

int getatm(char *fname, int maxnatoms, double **apos, char **atype, double *charge)
{
FILE *f;
int natoms;
	f = fopen(fname, "r");
        if(f==NULL) {
		fprintf(stderr,"can't open file %s\n",fname);
		return -1;
		}
	natoms = fgetatm(f, maxnatoms, apos, atype, charge);
	fclose(f);
	return natoms;
	}



int fgetatm(FILE *f, int maxnatoms, double **apos, char **atype, double *charge)
/*
 * Read a Brookhaven file putting the coords into apos and the char
 *  string names into atype. Return number of atoms.  M Pique
 * Stop reading on End-of-file, or 'END' line, or line with just a ".".
 * Here is some sample input: 
  1 0x7042  Ca  0.0335223  5.4441102  0.0069856   20    0   sp    Ca    0x3
  2 0x7042  Ca -2.3324949  8.1159780  1.7203546   20    0   sp    Ca    0x3
  3 0x7042   P  2.0533851  3.0058572  1.7203546   15    0  sp3     P    0x3
  END
 */
{
  char linebuf[200];
  register int n; /* number of atoms */

  char temp[7];
  int atnumb;

  n=0;
  while(fgets(linebuf, sizeof linebuf, f)!=NULL && n < maxnatoms &&
	0!=strncmp("END", linebuf,3) &&
	0!=strncmp("end", linebuf,3))
    {
      sscanf(linebuf,"%d %s %3s %lf %lf %lf", &atnumb, temp, atype[n],
	     &apos[n][0], &apos[n][1], &apos[n][2]);
#ifdef DEBUG
      printf("atom %d type %s %s at %lf %lf %lf\n",
	     n, atype[n], apos[n][0],  apos[n][1], apos[n][2]);
#endif
      charge[n] = 0.00;
      n++;
    }
  return n;
}

/* makes the atomtype indices */

int
makeatomtype(char **atype, int natoms, int *atomtype)
{

    /*
     * Generate into atomtype the integer 0..7 atom type Returns number of
     * atoms	
     *
     */

#include <ctype.h>
    register int i;
    char      a;
    char      b;
    int       anum;

    for (i = 0; i < natoms; i++) {
	a = atype[i][0];
	if (islower(a))
	    a = toupper(a);
	b = atype[i][1];
	if (islower(b))
	    b = toupper(b);
	switch (a) {
	  case 'H':
	    anum = 7;
	    break;
	  case 'C':
	    anum = 2;
	    if (b == 'A')
	      anum = 9;   /* calcium */
	    break;
	  case 'O':
	    anum = 4;
	    break;
	  case 'N':
	    anum = 3;
	    break;
	  case 'S':
	    anum = 6;
	    break;
	  case 'P':
	    anum = 8;
	    break;
	  case 'F':
	    anum = 10;
	    break;
	  default:
	    anum = 5;
	    break;
	}
	atomtype[i] = anum;
    }
    return natoms;
}

int makebonds(register double **apos, char **atype, int natoms, 
	      int maxnbonds, double threshold, int **bonds)
{
    /*
     * Generate into bonds the indices into apos of neighbors. Returns number
     * of bonds. 
     *
     * This version is a simple-minded distance check with special case code to
     * prevent hydrogen over-connectivity. Mike Pique 
     */

    register int i, j;
    register int nbonds;
    register int nbonds_this_atom;
    register double dx, dy, dxy, dz;

    nbonds = 0;
    for (i = natoms - 1; i > 0; i--) {
	nbonds_this_atom = 0;
	for (j = i - 1; j > 0 && nbonds < maxnbonds; j--) {
	    /*
	     * The outer loop index 'i' is AFTER the inner loop 'j': 'i'
	     * leads 'j' in the list: since hydrogens traditionally follow
	     * the heavy atom they're bonded to, this makes it easy to quit
	     * bonding to hydrogens after one bond is made by breaking out of
	     * the 'j' loop when 'i' is a hydrogen and we make a bond to it.
	     * Working backwards like this makes it easy to find the heavy
	     * atom that came 'just before' the Hydrogen. mp 
	     */

	    /* never bond hydrogens to each other... */
	    if (atype[i][0] == 'H' && atype[j][0] == 'H')
		continue;

	    dx = apos[i][0] - apos[j][0];
	    if (dx < 0)
		dx = -dx;
	    if (dx > threshold)
		continue;
	    dy = apos[i][1] - apos[j][1];
	    if (dy < 0)
		dy = -dy;
	    if (dy > threshold)
		continue;
/*	    if ((dxy = hypot(dx, dy)) > threshold) */
	    if ((dxy = sqrt((dx * dx) + (dy *dy))) > threshold)
		continue;

	    dz = apos[i][2] - apos[j][2];
	    if (dz < 0)
		dz = -dz;
	    if (dz > threshold)
		continue;
	    if (hypot(dxy, dz) > threshold)
		continue;

	    bonds[nbonds][0] = j;
	    bonds[nbonds][1] = i;
#ifdef DEBUG
	    printf("bond[%d] %d to %d\n", nbonds, i, j);
#endif
	    nbonds++;
	    nbonds_this_atom++;
	    if (atype[i][0] == 'H')
		break;		/* only one bond per hydrogen */
	    if (atype[i][0] == 'h')
		break;		/* only one bond per hydrogen */
	}
/*
	if (nbonds_this_atom == 0 && nbonds < maxnbonds) {
	    bonds[nbonds][0] = bonds[nbonds][1] = i;
	    nbonds++;
	}
*/
    }

    return nbonds;
}

int computebonds(register double **apos, char **atype, int natoms, double threshold)
{
    /*
     * Generate into bonds the indices into apos of neighbors. Returns number
     * of bonds. 
     *
     * This version is a simple-minded distance check with special case code to
     * prevent hydrogen over-connectivity. Mike Pique 
     */

    register int i, j;
    register int nbonds;
    register int nbonds_this_atom;
    register double dx, dy, dxy, dz;

    nbonds = 0;
    for (i = natoms - 1; i > 0; i--) {
	nbonds_this_atom = 0;
	for (j = i - 1; j > 0 ; j--) {
	    /*
	     * The outer loop index 'i' is AFTER the inner loop 'j': 'i'
	     * leads 'j' in the list: since hydrogens traditionally follow
	     * the heavy atom they're bonded to, this makes it easy to quit
	     * bonding to hydrogens after one bond is made by breaking out of
	     * the 'j' loop when 'i' is a hydrogen and we make a bond to it.
	     * Working backwards like this makes it easy to find the heavy
	     * atom that came 'just before' the Hydrogen. mp 
	     */

	    /* never bond hydrogens to each other... */
	    if (atype[i][0] == 'H' && atype[j][0] == 'H')
		continue;

	    dx = apos[i][0] - apos[j][0];
	    if (dx < 0)
		dx = -dx;
	    if (dx > threshold)
		continue;
	    dy = apos[i][1] - apos[j][1];
	    if (dy < 0)
		dy = -dy;
	    if (dy > threshold)
		continue;
	    if ((dxy = hypot(dx, dy)) > threshold)
		continue;

	    dz = apos[i][2] - apos[j][2];
	    if (dz < 0)
		dz = -dz;
	    if (dz > threshold)
		continue;
	    if (hypot(dxy, dz) > threshold)
		continue;

	    nbonds++;
	    nbonds_this_atom++;
	    if (atype[i][0] == 'H')
		break;		/* only one bond per hydrogen */
	    if (atype[i][0] == 'h')
		break;		/* only one bond per hydrogen */
	}
    }

    printf("found %d bonds\n", nbonds);
    return nbonds;
}


/* The main module begins */

main(int argc,char **argv)
{
  
  double fudge;
  double center_position[3];

  double **atoms = NULL;
  double **patoms = NULL;
  double *charges = NULL;
  int    *atom_types = NULL;
  int    **bonds = NULL;
  int    nbonds = 0;

  char **atype = NULL;
  char **amino = NULL;

  double max_scale;

  int do_pdb = 1;
  int do_sky = 0;
  int do_ground = 0;
  int do_check = 0;
  int do_cpk = 0;
  int do_cov = 0;
  int do_vdw = 1;
  int do_ball = 0;
  int do_hybrid = 0;
  int do_area = 0;
  int do_plain = 0;

  double rad_factor = 1.0;
  double threshold = 2.2;

  int natoms;
  int rotation_axis, rotation_frames, rotation_increment;
  int myerr = 0;

  int obj_only = 0;
  char inputfname[128],outputfname[128];
  
  if (0 == parse_inputs(argc,argv,inputfname,outputfname,&fudge,
			&rotation_axis,&rotation_frames,
			&rotation_increment,&rad_factor, 
			&obj_only, &do_pdb, &do_sky, &do_ground, &do_check,
			&do_vdw, &do_cov, &do_cpk, &do_ball, &threshold,
			&do_hybrid, &do_area, &do_plain)) 
    {
      return(ERR_PARSE_ARGS);
    }

  printf("Scanning  atom file <%s>... ", inputfname);
  natoms = getpdb_atnum(inputfname);

  if (natoms <= 0)
    {
      printf("Couldn't read any atoms from <%s>! Aborting.\n", inputfname);
      myerr = ERR_NO_ATOMS;
      exit(myerr);
    }

  printf("Got <%d> atoms. \n", natoms);

  if ((atoms = dmatrix(natoms, 3)) != NULL) {
    printf("dmatrix \n");
      if ((patoms = dmatrix(natoms, 3)) != NULL) {
	printf("patoms \n");
	  if ((charges = dvector(natoms)) != NULL) {
	    printf("charges \n");
	      if ((atype = cmatrix(natoms, 4)) != NULL) {
		printf("atype \n");
		  if ((amino = cmatrix(natoms, 4)) != NULL) {
		    printf("amino \n");
		      if ((atom_types = ivector(natoms)) != NULL) {
			  if (do_pdb)
			    getpdb(inputfname,natoms,atoms,atype,charges,amino);
			  else
			    getatm(inputfname,natoms,atoms,atype,charges);
  
			  makeatomtype(atype,natoms,atom_types);
			  max_scale = 1;
			  
			  printf("Re-orienting and positioning protein. \n");

			  rotate_protein(natoms,atoms,patoms,xrot, yrot, zrot); 
			  center_protein(natoms,patoms,center_position);
			  
			  flip_zaxis(natoms,patoms); /* for left handed coordinate system */
			  
			  if (do_ball)
			    {
			      printf("Computing bonds: ");
			      if ((nbonds = computebonds(patoms, atype, natoms, threshold)) <= 0)
				{
				  free_ivector(atom_types, natoms);
				  free_cmatrix(amino, natoms, 4);
				  free_cmatrix(atype, natoms, 4);
				  free_dvector(charges, natoms);
				  free_dmatrix(patoms, natoms, 3);
				  free_dmatrix(atoms, natoms, 3);
				  myerr = ERR_NO_BONDS;
				  printf("No bonds found? Errorcode is: %d\n", ERR_NO_BONDS);
				  exit(myerr);
				}
			      if ((bonds = imatrix(nbonds,2))!=NULL)
				{
				  makebonds(patoms, atype, natoms, nbonds, threshold, bonds);
				  rad_factor *= .3;
				  write_output(outputfname, natoms, patoms, atom_types, 
					       bonds, nbonds, rad_factor, obj_only, 
					       do_sky, do_ground, do_check, do_cpk, 
					       do_cov, do_vdw, do_ball, do_hybrid, do_area, do_plain);
				  free_imatrix(bonds,nbonds,2);
				}
			      else
				{
				  free_ivector(atom_types, natoms);
				  free_cmatrix(amino, natoms, 4);
				  free_cmatrix(atype, natoms, 4);
				  free_dvector(charges, natoms);
				  free_dmatrix(patoms, natoms, 3);
				  free_dmatrix(atoms, natoms, 3);
				  printf("can't allocate memory for bonds??\n");
				  myerr = ERR_CANT_ALLOC_BONDS;
				  exit(myerr);
				}
			    }
			  else
			    write_output(outputfname, natoms, patoms, atom_types, 
					 bonds, nbonds, rad_factor, obj_only, 
					 do_sky, do_ground, do_check, do_cpk, 
					 do_cov, do_vdw, do_ball, do_hybrid, do_area, do_plain);
			  
			  free_ivector(atom_types, natoms);
			}
		      free_cmatrix(amino, natoms, 4);
		    }
		  free_cmatrix(atype, natoms, 4);
		}
	      free_dvector(charges, natoms);
	    }
	  free_dmatrix(patoms, natoms, 3);
	}
      free_dmatrix(atoms, natoms, 3);
    }
  exit(0);
}

/******************************* End of pdb2pov.c ****************************/

