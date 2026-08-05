/*
 * pdb2pov.h -- shared constants, types and prototypes for pdb2pov.
 *
 * Replaces the 1993 pdb2light.h, which carried a pile of constants that were
 * never referenced plus two conflicting definitions of X/Y/Z.
 *
 * Copyright (c) 1993-2026 Eric G. Suchanek, Ph.D.  All rights reserved.
 * Subject to the GNU License.
 */

#ifndef PDB2POV_H
#define PDB2POV_H

#include <stdio.h>

/* --- version ----------------------------------------------------------- */

#define PDB2POV_VERSION "2.1"

/*
 * POV-Ray language level the emitted scenes are written against.  3.7 is the
 * stable release carried by the distributions; it requires every #declare to
 * end in a semicolon and every scene to set assumed_gamma, both of which the
 * writer now does.
 */
#define POV_VERSION "3.7"

/* --- coordinate indices ------------------------------------------------- */

#define X 0
#define Y 1
#define Z 2

/* --- atom type codes ---------------------------------------------------- */

/*
 * Indices produced by make_atom_types().  The values are historical and are
 * not contiguous; ANY_TYPE is the catch-all for elements with no dedicated
 * texture in atoms2.inc.
 */
#define H_TYPE   7
#define C_TYPE   2
#define O_TYPE   4
#define N_TYPE   3
#define S_TYPE   6
#define P_TYPE   8
#define CA_TYPE  9
#define FE_TYPE 10
#define ANY_TYPE 5

/* --- geometry ----------------------------------------------------------- */

/*
 * Largest atomic radius in angstroms, used to pad the bounding box and the
 * enclosing sphere.  Deliberately a generous overestimate: the padding is
 * applied to unscaled coordinates, so it stays conservative once ATM_SCL
 * shrinks the spheres in ball-and-stick mode.
 */
#define MAX_RAD 2.5

/* Radius of the cylinders drawn between bonded atoms, in angstroms. */
#define BOND_RAD 0.17

/* Bounding sphere is grown by this fraction so it clears the outermost atom. */
#define SPHERE_FUDGE 0.02

/* --- buffer sizes ------------------------------------------------------- */

#define PATH_MAX_LEN 1024 /* input/output paths, extension included  */
#define LINE_MAX_LEN 256  /* one record from a PDB or .atm file      */
#define ATOM_NAME_LEN 5   /* PDB atom name field, plus terminator    */
#define RES_NAME_LEN 5    /* PDB residue name field, plus terminator */
#define ELEM_NAME_LEN 3   /* PDB element symbol, plus terminator     */
#define CHAIN_FILTER_LEN 32 /* accepted chain IDs, as a string       */

/* --- rendering options -------------------------------------------------- */

/*
 * Which set of atomic radii the scene includes.
 *
 * No command-line flag selects RADII_CPK, and none did in 1993 either -- the
 * original tracked a do_cpk flag that nothing ever set.  It is kept because
 * atoms_cpk.inc still ships, but note that its radii are identical to
 * atoms_vdw.inc's, so the option would be a no-op if it were wired up.
 */
typedef enum {
    RADII_VDW = 0,      /* atoms_vdw.inc      */
    RADII_COVALENT,     /* atoms_covalent.inc */
    RADII_CPK           /* atoms_cpk.inc      */
} RadiiSet;

/* What, if anything, surrounds the molecule. */
typedef enum {
    BACKDROP_NONE = 0,  /* neither sky nor ground   */
    BACKDROP_SKY,       /* gradient sky with clouds */
    BACKDROP_PLAIN      /* flat white sky           */
} Backdrop;

typedef enum {
    GROUND_NONE = 0,
    GROUND_PLAIN,       /* bumped RichBlue plane */
    GROUND_CHECKER      /* the classic checkerboard */
} GroundStyle;

/* Everything parse_options() collects from argv. */
typedef struct {
    char input_path[PATH_MAX_LEN];  /* with .pdb or .atm appended */
    char output_stem[PATH_MAX_LEN]; /* without extension          */

    double xrot, yrot, zrot;        /* absolute rotations, degrees */
    double radii_scale;             /* ATM_SCL before ball-mode shrink */
    double bond_threshold;          /* bond cutoff, angstroms */

    RadiiSet radii;
    Backdrop backdrop;
    GroundStyle ground;

    int atm_format;   /* read .atm instead of .pdb        */
    int object_only;  /* emit a bare .inc, no camera/light */
    int ball_stick;   /* draw bond cylinders              */
    int glass_atoms;  /* overlay the glass-atom merge     */
    int area_light;   /* soft light instead of a point    */

    /*
     * Parser behaviour.  The defaults changed in 2.1; see read_pdb().
     */
    int legacy_elements; /* guess elements from atom names, pre-2.1 style */
    int keep_altlocs;    /* keep every alternate conformation             */
    char chain_filter[CHAIN_FILTER_LEN]; /* accepted chain IDs; "" = all  */
} Options;

/*
 * What a parse discarded, so it can be reported rather than swallowed.  The
 * 1993 code dropped atoms silently in three separate places.
 */
#define MAX_REPORTED_SYMBOLS 12

typedef struct {
    int accepted;
    int skipped_altloc;
    int skipped_chain;
    int skipped_malformed;
    int freeform_fallback; /* records whose columns would not parse */
    int generic;           /* atoms with no dedicated texture       */
    int no_element_column; /* records with columns 77-78 blank      */
    char generic_symbols[MAX_REPORTED_SYMBOLS][ELEM_NAME_LEN];
    int n_generic_symbols;
} ParseStats;

/* A parsed structure: coordinates plus the per-atom fields we keep. */
typedef struct {
    int natoms;
    double **pos;     /* [natoms][3] */
    char **atom_name; /* [natoms][ATOM_NAME_LEN] */
    char **res_name;  /* [natoms][RES_NAME_LEN]  */
    char **element;   /* [natoms][ELEM_NAME_LEN], upper case, "" if unknown */
    double *charge;   /* [natoms] */
    int *type;        /* [natoms], one of *_TYPE */
} Structure;

/* An axis-aligned bounding box, already padded by MAX_RAD. */
typedef struct {
    double xmin, xmax;
    double ymin, ymax;
    double zmin, zmax;
} Extents;

/* --- prototypes --------------------------------------------------------- */

/* pdb_io.c-equivalent routines, all still in pdb2pov.c */
int count_pdb_atoms(const char *path);
int read_pdb(const char *path, int max_atoms, const Options *opt, Structure *s,
             ParseStats *st);
int read_atm(const char *path, int max_atoms, Structure *s, ParseStats *st);
void make_atom_types(const Options *opt, Structure *s, ParseStats *st);
void report_parse(const ParseStats *st, int legacy);

/* geometry */
void make_rotmat(double rotmat[3][3], double xdeg, double ydeg, double zdeg);
void matrix_times_vector(const double m[3][3], double *x, double *y, double *z);
void rotate_structure(Structure *s, double xdeg, double ydeg, double zdeg);
void center_structure(Structure *s, double center[3]);
void flip_zaxis(Structure *s);
void compute_extents(const Structure *s, Extents *e);
double compute_sphere(const Structure *s);

/* bonds */
int find_bonds(const Structure *s, double threshold, int max_bonds,
               int **bonds);

/* output */
int write_output(const Options *opt, const Structure *s, int **bonds,
                 int nbonds);

#endif /* PDB2POV_H */
