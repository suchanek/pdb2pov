/*
 * pdb2pov.c -- convert Brookhaven PDB atomic structure files into POV-Ray
 * scene files.
 *
 * Author: Eric G. Suchanek, Ph.D.
 * Copyright (c) 1993-2026 Eric G. Suchanek, Ph.D.  All rights reserved.
 * Subject to the GNU License.
 *
 * The emitted scenes require the bundled radius and texture includes
 * (atoms_vdw.inc, atoms2.inc and friends) to be on POV-Ray's library path.
 *
 * Modernisation notes (2026).  The 1993 sources were K&R C with AmigaOS
 * branches, and needed a set of Makefile force-includes to build at all on a
 * 64-bit host.  What changed:
 *
 *   - K&R definitions became prototypes, so the compiler can check calls.
 *     This removes the pointer-truncation class of bug that the old
 *     pdb2pov_protos.h force-include existed to work around.
 *   - The AmigaOS branches, asyncio.c and m68040.c are gone.  Git history
 *     preserves them.
 *   - Fixed-size string handling moved to snprintf.  The old date formatting
 *     overran its buffer -- which is why the Makefile had to disable
 *     _FORTIFY_SOURCE -- and printed a two-digit year that read "126" in 2026.
 *   - The bond search skipped atom 0 entirely (`j > 0` where `j >= 0` was
 *     meant), so the first atom in every file was never bonded.
 *   - The scene writer emits POV-Ray 3.7 rather than 2.x: semicolons after
 *     float declarations, commas in vector/scalar pairs, light sources no
 *     longer wrapped in object{}, refraction moved from finish to interior,
 *     and colour_map ramps rewritten as 3.x control points.
 *
 * Scene geometry is deliberately unchanged -- see DEG_PER_RAD_* below.
 */

#include <ctype.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "pdb2pov.h"
#include "pdb2pov_errors.h"
#include "pdb2pov_usage.h"
#include "util.h"

/*
 * The 1993 code converted degrees to radians with two hand-rolled constants:
 * 57.29 in the rotation matrix, and 57.298 when placing the camera and light.
 * Both are truncations of 180/pi (57.29577951...).
 *
 * They are kept verbatim rather than replaced with an exact conversion.  The
 * scene geometry these produce is quoted downstream -- quiltwright's
 * docs/pdb2pov.md pins crambin's camera distance at 40.075 and its enclosing
 * sphere at 18.759, both derived from these exact values.  Correcting the
 * constants would silently move every camera in every existing scene for a
 * change of well under a tenth of a degree.
 */
#define DEG_PER_RAD_ROT 57.29  /* make_rotmat()                     */
#define DEG_PER_RAD_CAM 57.298 /* write_camera(), write_light()     */

#define CAMERA_HALF_ANGLE_DEG 22.0
#define LIGHT_HALF_ANGLE_DEG 20.0

/* ------------------------------------------------------------------ */
/* Option parsing                                                     */
/* ------------------------------------------------------------------ */

static void options_init(Options *opt)
{
    memset(opt, 0, sizeof *opt);

    /*
     * Historical defaults: van der Waals radii, no backdrop, no ground,
     * space-filling (no bond cylinders), 2.2 angstrom bond cutoff.
     */
    opt->radii = RADII_VDW;
    opt->backdrop = BACKDROP_NONE;
    opt->ground = GROUND_NONE;
    opt->radii_scale = 1.0;
    opt->bond_threshold = 2.2;
}

/*
 * Read a double from the argument following argv[i], advancing i.  Returns 0
 * and reports if the flag was last on the line or the value does not parse --
 * the old code ran off the end of argv and scanned whatever followed.
 */
static int take_double(int argc, char **argv, int *i, double *out,
                       const char *flag)
{
    char *end;
    double v;

    if (*i + 1 >= argc) {
        fprintf(stderr, "pdb2pov: %s requires a value\n", flag);
        return 0;
    }

    v = strtod(argv[++(*i)], &end);
    if (end == argv[*i] || *end != '\0') {
        fprintf(stderr, "pdb2pov: %s: '%s' is not a number\n", flag,
                argv[*i]);
        return 0;
    }

    *out = v;
    return 1;
}

/*
 * Build "stem.ext" in dst, reporting rather than truncating silently.  The
 * original used strcpy/strcat into a 128-byte buffer with no length check.
 */
static int join_extension(char *dst, size_t dstsz, const char *stem,
                          const char *ext)
{
    int n = snprintf(dst, dstsz, "%s%s", stem, ext);

    if (n < 0 || (size_t)n >= dstsz) {
        fprintf(stderr, "pdb2pov: path '%s%s' is too long\n", stem, ext);
        return 0;
    }
    return 1;
}

/*
 * Read the argument following argv[i] as a string, advancing i.
 */
static int take_string(int argc, char **argv, int *i, char *out, size_t outsz,
                       const char *flag)
{
    int n;

    if (*i + 1 >= argc) {
        fprintf(stderr, "pdb2pov: %s requires a value\n", flag);
        return 0;
    }

    n = snprintf(out, outsz, "%s", argv[++(*i)]);
    if (n < 0 || (size_t)n >= outsz) {
        fprintf(stderr, "pdb2pov: %s: value is too long\n", flag);
        return 0;
    }
    return 1;
}

/*
 * Long options carry the parser behaviour added in 2.1.  They are spelled out
 * rather than given letters because the single-character space is full: -c is
 * covalent radii and -a is the area light, so --chain and --keep-altlocs have
 * nowhere obvious to go.
 */
static int parse_long_option(int argc, char **argv, int *i, Options *opt)
{
    const char *arg = argv[*i];

    if (strcmp(arg, "--legacy-elements") == 0) {
        opt->legacy_elements = 1;
        return 1;
    }
    if (strcmp(arg, "--keep-altlocs") == 0) {
        opt->keep_altlocs = 1;
        return 1;
    }
    if (strcmp(arg, "--chain") == 0) {
        return take_string(argc, argv, i, opt->chain_filter,
                           sizeof opt->chain_filter, "--chain");
    }

    fprintf(stderr, "pdb2pov: unrecognised option '%s'\n", arg);
    return 0;
}

static int parse_options(int argc, char **argv, Options *opt)
{
    int i;

    options_init(opt);

    if (argc < 3) {
        fprintf(stderr, PDB2POV_USAGE, argv[0], argv[0]);
        return 0;
    }

    for (i = 3; i < argc && argv[i][0] == '-'; i++) {
        if (strncmp(argv[i], "--", 2) == 0) {
            if (!parse_long_option(argc, argv, &i, opt))
                return 0;
            continue;
        }

        switch (argv[i][1]) {
        case 'v':
        case 'V':
            opt->radii = RADII_VDW;
            break;
        case 'c':
        case 'C':
            opt->radii = RADII_COVALENT;
            break;
        case 'b':
        case 'B':
            opt->ball_stick = 1;
            break;
        case 'q':
        case 'Q':
            opt->ball_stick = 1;
            opt->glass_atoms = 1;
            break;
        case 'd':
        case 'D':
            if (!take_double(argc, argv, &i, &opt->bond_threshold, "-d"))
                return 0;
            break;
        case 'r':
        case 'R':
            if (!take_double(argc, argv, &i, &opt->radii_scale, "-r"))
                return 0;
            if (opt->radii_scale <= 0.0)
                opt->radii_scale = 1.0;
            break;
        case 'x':
        case 'X':
            if (!take_double(argc, argv, &i, &opt->xrot, "-x"))
                return 0;
            break;
        case 'y':
        case 'Y':
            if (!take_double(argc, argv, &i, &opt->yrot, "-y"))
                return 0;
            break;
        case 'z':
        case 'Z':
            if (!take_double(argc, argv, &i, &opt->zrot, "-z"))
                return 0;
            break;
        case 'o':
        case 'O':
            opt->object_only = 1;
            break;
        case 't':
        case 'T':
            opt->atm_format = 1;
            break;
        case 'a':
        case 'A':
            opt->area_light = 1;
            break;
        case 'p':
        case 'P':
            opt->backdrop = BACKDROP_PLAIN;
            opt->ground = GROUND_NONE;
            break;
        case 's':
        case 'S':
            opt->backdrop = BACKDROP_SKY;
            break;
        case 'g':
        case 'G':
            opt->ground = GROUND_PLAIN;
            break;
        case 'h':
        case 'H':
            opt->ground = GROUND_CHECKER;
            break;
        default:
            fprintf(stderr, "pdb2pov: unrecognised option '%s'\n", argv[i]);
            fprintf(stderr, PDB2POV_USAGE, argv[0], argv[0]);
            return 0;
        }
    }

    if (i < argc) {
        fprintf(stderr, "pdb2pov: unexpected argument '%s'\n", argv[i]);
        return 0;
    }

    if (!join_extension(opt->input_path, sizeof opt->input_path, argv[1],
                        opt->atm_format ? ".atm" : ".pdb"))
        return 0;

    if (!join_extension(opt->output_stem, sizeof opt->output_stem, argv[2],
                        ""))
        return 0;

    return 1;
}

/* ------------------------------------------------------------------ */
/* PDB and .atm input                                                 */
/* ------------------------------------------------------------------ */

/* True for the record types that carry coordinates. */
static int is_coordinate_record(const char *line)
{
    return strncmp(line, "ATOM", 4) == 0 || strncmp(line, "atom", 4) == 0 ||
           strncmp(line, "HETATM", 6) == 0 || strncmp(line, "hetatm", 6) == 0;
}

/*
 * Copy a fixed PDB column range, trim surrounding blanks, and upper-case it.
 * Returns the trimmed length.  PDB fields are blank-padded and some are
 * right-justified (the element symbol) while others are not, so trimming both
 * ends is the only thing that works across fields.
 */
static size_t field_text(const char *line, size_t len, size_t start,
                         size_t width, char *out, size_t outsz, int upper)
{
    size_t i, n = 0;

    out[0] = '\0';
    if (start >= len)
        return 0;
    if (start + width > len)
        width = len - start;

    for (i = 0; i < width && n + 1 < outsz; i++) {
        unsigned char c = (unsigned char)line[start + i];

        if (c == '\n' || c == '\r')
            break;
        out[n++] = (char)(upper ? toupper(c) : c);
    }
    out[n] = '\0';

    /* trim trailing blanks */
    while (n > 0 && (out[n - 1] == ' ' || out[n - 1] == '\t'))
        out[--n] = '\0';

    /* trim leading blanks */
    i = 0;
    while (out[i] == ' ' || out[i] == '\t')
        i++;
    if (i > 0) {
        memmove(out, out + i, n - i + 1);
        n -= i;
    }

    return n;
}

/* Parse a fixed column range as a double.  Returns 0 if it is not a number. */
static int field_double(const char *line, size_t len, size_t start,
                        size_t width, double *out)
{
    char buf[32];
    char *end;
    double v;

    if (field_text(line, len, start, width, buf, sizeof buf, 0) == 0)
        return 0;

    v = strtod(buf, &end);
    if (end == buf || *end != '\0')
        return 0;

    *out = v;
    return 1;
}

/*
 * Map a PDB element symbol to one of the *_TYPE codes.  Returns ANY_TYPE for
 * anything without a dedicated texture in atoms2.inc; such atoms render with
 * the generic Atom_X rather than being dropped.
 */
static int element_to_type(const char *el)
{
    if (strcmp(el, "H") == 0 || strcmp(el, "D") == 0)
        return H_TYPE; /* deuterium renders as hydrogen */
    if (strcmp(el, "C") == 0)
        return C_TYPE;
    if (strcmp(el, "N") == 0)
        return N_TYPE;
    if (strcmp(el, "O") == 0)
        return O_TYPE;
    if (strcmp(el, "S") == 0)
        return S_TYPE;
    if (strcmp(el, "P") == 0)
        return P_TYPE;
    if (strcmp(el, "CA") == 0)
        return CA_TYPE;
    if (strcmp(el, "FE") == 0)
        return FE_TYPE;
    return ANY_TYPE;
}

/*
 * The pre-2.1 element guess: look at the first one or two characters of the
 * atom name.  Kept for --legacy-elements, and used as a fallback for files
 * old enough to predate the element column.
 *
 * It is wrong for any two-letter element whose first letter collides with a
 * one-letter element -- sodium reads as nitrogen, chlorine as carbon,
 * fluorine as iron -- and it cannot tell a protein alpha carbon from calcium
 * without the caller's help.
 */
static int guess_type_from_name(const char *name)
{
    int a = toupper((unsigned char)name[0]);
    int b = toupper((unsigned char)name[1]);

    switch (a) {
    case 'H': return H_TYPE;
    case 'C': return (b == 'A') ? CA_TYPE : C_TYPE;
    case 'O': return O_TYPE;
    case 'N': return N_TYPE;
    case 'S': return S_TYPE;
    case 'P': return P_TYPE;
    case 'F': return FE_TYPE;
    default:  return ANY_TYPE;
    }
}

/* True for the records that terminate a coordinate section. */
static int is_terminator(const char *line)
{
    return strcmp(line, ".") == 0 || strncmp(line, "END", 3) == 0 ||
           strncmp(line, "end", 3) == 0;
}

/*
 * Count coordinate records so the caller can size its arrays.  Returns the
 * count, or a negative ERR_* code.
 */
int count_pdb_atoms(const char *path)
{
    char line[LINE_MAX_LEN];
    FILE *f;
    int n = 0;

    f = fopen(path, "r");
    if (f == NULL) {
        fprintf(stderr, "pdb2pov: can't open '%s' for reading\n", path);
        return ERR_CANT_READ_INPUT;
    }

    while (fgets(line, sizeof line, f) != NULL && !is_terminator(line)) {
        if (is_coordinate_record(line))
            n++;
    }

    fclose(f);
    return n;
}

/*
 * Read a Brookhaven file into s.  Originally by M. Pique; the column offsets
 * are unchanged.
 *
 * Sample input:
 *   ATOM     21  N   CYS     3      13.507  11.238   8.398
 *   HETATM   99  FE  XXX     4      33.265  10.443  10.307
 *   HETATM       CU+2CU    152     -20.577   2.601 -14.587       2.000
 *
 * The trailing charge column is not genuine PDB and is optional.
 *
 * The offsets are reached with an explicit length check now.  The old code
 * indexed linebuf[70] unconditionally once a record matched, reading past the
 * terminator on any line shorter than that -- which is most of them.
 */
int read_pdb(const char *path, int max_atoms, const Options *opt, Structure *s,
             ParseStats *st)
{
    char line[LINE_MAX_LEN];
    FILE *f;
    int n = 0;

    f = fopen(path, "r");
    if (f == NULL) {
        fprintf(stderr, "pdb2pov: can't open '%s' for reading\n", path);
        return ERR_CANT_READ_INPUT;
    }

    while (fgets(line, sizeof line, f) != NULL && n < max_atoms &&
           !is_terminator(line)) {
        size_t len = strlen(line);
        double charge;
        char alt_loc, chain_id;

        if (!is_coordinate_record(line))
            continue;

        /* Coordinates start at column 31; anything shorter is not a record. */
        if (len <= 30) {
            st->skipped_malformed++;
            continue;
        }

        /*
         * Alternate conformations.  A residue modelled two ways contributes
         * both sets of atoms, which renders as overlapping spheres plus
         * spurious bonds between the A and B conformers.  Keep the blank and
         * 'A' indicators only, unless asked otherwise.
         */
        alt_loc = (len > 16) ? line[16] : ' ';
        if (!opt->keep_altlocs && alt_loc != ' ' && alt_loc != 'A' &&
            alt_loc != 'a') {
            st->skipped_altloc++;
            continue;
        }

        chain_id = (len > 21) ? line[21] : ' ';
        if (opt->chain_filter[0] != '\0' &&
            strchr(opt->chain_filter, chain_id) == NULL) {
            st->skipped_chain++;
            continue;
        }

        /* Atom name, columns 13-16. */
        if (field_text(line, len, 12, 4, s->atom_name[n], ATOM_NAME_LEN, 0) == 0) {
            st->skipped_malformed++;
            continue;
        }

        /* Residue name, columns 18-20.  Absence is tolerated. */
        field_text(line, len, 17, 3, s->res_name[n], RES_NAME_LEN, 0);

        /*
         * Coordinates, columns 31-38 / 39-46 / 47-54.  Fixed columns are what
         * the format actually specifies.  If a record does not honour them --
         * hand-edited files exist -- fall back to the free-form scan the
         * pre-2.1 parser used, so nothing that converted before stops
         * converting now.
         */
        if (!field_double(line, len, 30, 8, &s->pos[n][X]) ||
            !field_double(line, len, 38, 8, &s->pos[n][Y]) ||
            !field_double(line, len, 46, 8, &s->pos[n][Z])) {
            if (sscanf(&line[30], "%lf %lf %lf", &s->pos[n][X], &s->pos[n][Y],
                       &s->pos[n][Z]) != 3) {
                st->skipped_malformed++;
                continue;
            }
            st->freeform_fallback++;
        }

        /* Element symbol, columns 77-78.  Blank on files that predate it. */
        field_text(line, len, 76, 2, s->element[n], ELEM_NAME_LEN, 1);
        if (s->element[n][0] == '\0')
            st->no_element_column++;

        /*
         * The alpha-carbon hack, needed only when guessing from the name: in
         * an ATOM record "CA" is C-alpha, not calcium.  Blanking the second
         * character makes the guess read it as carbon.  With a real element
         * column the distinction is free, so this is confined to the legacy
         * path.
         */
        if ((opt->legacy_elements || s->element[n][0] == '\0') &&
            strncmp(line, "ATOM", 4) == 0 &&
            (strncmp(s->atom_name[n], "CA", 2) == 0 ||
             strncmp(s->atom_name[n], "ca", 2) == 0))
            s->atom_name[n][1] = ' ';

        /* Non-standard trailing charge column, if present. */
        if (len > 70 && sscanf(&line[70], "%lf", &charge) == 1)
            s->charge[n] = charge;
        else
            s->charge[n] = 0.0;

        n++;
    }

    fclose(f);
    s->natoms = n;
    st->accepted = n;
    return n;
}

/*
 * Read the alternative .atm format:
 *   1 0x7042  Ca  0.0335223  5.4441102  0.0069856   20    0   sp    Ca    0x3
 */
int read_atm(const char *path, int max_atoms, Structure *s, ParseStats *st)
{
    char line[LINE_MAX_LEN];
    char flags[32];
    FILE *f;
    int n = 0;
    int atom_number;

    f = fopen(path, "r");
    if (f == NULL) {
        fprintf(stderr, "pdb2pov: can't open '%s' for reading\n", path);
        return ERR_CANT_READ_INPUT;
    }

    while (fgets(line, sizeof line, f) != NULL && n < max_atoms &&
           !is_terminator(line)) {
        if (sscanf(line, "%d %31s %3s %lf %lf %lf", &atom_number, flags,
                   s->atom_name[n], &s->pos[n][X], &s->pos[n][Y],
                   &s->pos[n][Z]) != 6) {
            st->skipped_malformed++;
            continue;
        }

        s->res_name[n][0] = '\0';
        /* The .atm format has no element column; types come from the name. */
        s->element[n][0] = '\0';
        s->charge[n] = 0.0;
        n++;
    }

    fclose(f);
    s->natoms = n;
    st->accepted = n;
    /* Not a defect for this format, so do not report it as one. */
    st->no_element_column = 0;
    return n;
}

/* Record an element symbol that fell through to the generic texture. */
static void note_generic_symbol(ParseStats *st, const char *el)
{
    int i;

    if (el[0] == '\0')
        el = "?";

    for (i = 0; i < st->n_generic_symbols; i++) {
        if (strcmp(st->generic_symbols[i], el) == 0)
            return;
    }

    if (st->n_generic_symbols < MAX_REPORTED_SYMBOLS) {
        snprintf(st->generic_symbols[st->n_generic_symbols], ELEM_NAME_LEN,
                 "%s", el);
        st->n_generic_symbols++;
    }
}

/*
 * Assign each atom one of the *_TYPE codes, which select a texture in
 * atoms2.inc.
 *
 * The element column drives this by default.  Before 2.1 the type was guessed
 * from the atom name, which mistyped every two-letter element sharing a first
 * letter with a one-letter one (NA as nitrogen, CL as carbon, F as iron) and
 * dropped the rest without a word.  --legacy-elements restores that.
 *
 * Anything with no dedicated texture is ANY_TYPE and renders as the generic
 * Atom_X, so an unrecognised element can no longer vanish from a scene.
 */
void make_atom_types(const Options *opt, Structure *s, ParseStats *st)
{
    int i;

    for (i = 0; i < s->natoms; i++) {
        int type;

        if (opt->legacy_elements || s->element[i][0] == '\0')
            type = guess_type_from_name(s->atom_name[i]);
        else
            type = element_to_type(s->element[i]);

        s->type[i] = type;

        if (type == ANY_TYPE) {
            st->generic++;
            /*
             * Report the element symbol when there is one; otherwise the atom
             * name is the only identifying thing the record carried.
             */
            note_generic_symbol(st, s->element[i][0] != '\0' && !opt->legacy_elements
                                        ? s->element[i]
                                        : s->atom_name[i]);
        }
    }
}

/* Say what the parse discarded or had to guess at. */
void report_parse(const ParseStats *st, int legacy)
{
    int i;

    if (st->skipped_altloc)
        printf("  %d atom(s) in alternate conformations skipped "
               "(--keep-altlocs to keep them)\n", st->skipped_altloc);
    if (st->skipped_chain)
        printf("  %d atom(s) outside the requested chain(s) skipped\n",
               st->skipped_chain);
    if (st->skipped_malformed)
        printf("  %d unparseable coordinate record(s) skipped\n",
               st->skipped_malformed);
    if (st->freeform_fallback)
        printf("  %d record(s) did not honour the PDB column layout; "
               "read free-form\n", st->freeform_fallback);
    if (st->no_element_column)
        printf("  %d record(s) have no element column; guessed from the atom "
               "name\n", st->no_element_column);

    if (st->generic) {
        printf("  %d atom(s) have no dedicated texture and %s (",
               st->generic,
               legacy ? "are dropped (--legacy-elements)"
                      : "render as Atom_X");
        for (i = 0; i < st->n_generic_symbols; i++)
            printf("%s%s", i ? ", " : "", st->generic_symbols[i]);
        if (st->n_generic_symbols == MAX_REPORTED_SYMBOLS)
            printf(", ...");
        printf(")\n");
    }
}

/* ------------------------------------------------------------------ */
/* Geometry                                                           */
/* ------------------------------------------------------------------ */

void make_rotmat(double rotmat[3][3], double xdeg, double ydeg, double zdeg)
{
    double g = xdeg / DEG_PER_RAD_ROT;
    double b = ydeg / DEG_PER_RAD_ROT;
    double t = zdeg / DEG_PER_RAD_ROT;

    double sing = sin(g), cosg = cos(g);
    double sinb = sin(b), cosb = cos(b);
    double sint = sin(t), cost = cos(t);

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

void matrix_times_vector(const double m[3][3], double *x, double *y, double *z)
{
    double xx = m[0][0] * *x + m[0][1] * *y + m[0][2] * *z;
    double yy = m[1][0] * *x + m[1][1] * *y + m[1][2] * *z;
    double zz = m[2][0] * *x + m[2][1] * *y + m[2][2] * *z;

    *x = xx;
    *y = yy;
    *z = zz;
}

/* Apply absolute rotations about each axis, in place. */
void rotate_structure(Structure *s, double xdeg, double ydeg, double zdeg)
{
    double matrix[3][3];
    int i;

    if (xdeg == 0.0 && ydeg == 0.0 && zdeg == 0.0)
        return;

    make_rotmat(matrix, xdeg, ydeg, zdeg);

    for (i = 0; i < s->natoms; i++)
        /*
         * The cast is required because C before C2x does not implicitly
         * convert double(*)[3] to const double(*)[3].
         */
        matrix_times_vector((const double (*)[3])matrix, &s->pos[i][X],
                            &s->pos[i][Y], &s->pos[i][Z]);
}

/* Translate to the geometric centre of coordinates, reporting the shift. */
void center_structure(Structure *s, double center[3])
{
    double cx = 0.0, cy = 0.0, cz = 0.0;
    int i;

    if (s->natoms == 0) {
        center[X] = center[Y] = center[Z] = 0.0;
        return;
    }

    for (i = 0; i < s->natoms; i++) {
        cx += s->pos[i][X];
        cy += s->pos[i][Y];
        cz += s->pos[i][Z];
    }

    cx /= s->natoms;
    cy /= s->natoms;
    cz /= s->natoms;

    for (i = 0; i < s->natoms; i++) {
        s->pos[i][X] -= cx;
        s->pos[i][Y] -= cy;
        s->pos[i][Z] -= cz;
    }

    center[X] = cx;
    center[Y] = cy;
    center[Z] = cz;
}

/* PDB coordinates are right-handed; POV-Ray's world is left-handed. */
void flip_zaxis(Structure *s)
{
    int i;

    for (i = 0; i < s->natoms; i++)
        s->pos[i][Z] = -s->pos[i][Z];
}

/* Axis-aligned bounds, padded by MAX_RAD to allow for the atom spheres. */
void compute_extents(const Structure *s, Extents *e)
{
    double xmin = 9999.9, ymin = 9999.9, zmin = 9999.9;
    double xmax = -9999.9, ymax = -9999.9, zmax = -9999.9;
    int i;

    for (i = 0; i < s->natoms; i++) {
        if (s->pos[i][X] < xmin) xmin = s->pos[i][X];
        if (s->pos[i][X] > xmax) xmax = s->pos[i][X];
        if (s->pos[i][Y] < ymin) ymin = s->pos[i][Y];
        if (s->pos[i][Y] > ymax) ymax = s->pos[i][Y];
        if (s->pos[i][Z] < zmin) zmin = s->pos[i][Z];
        if (s->pos[i][Z] > zmax) zmax = s->pos[i][Z];
    }

    e->xmin = xmin - MAX_RAD;
    e->xmax = xmax + MAX_RAD;
    e->ymin = ymin - MAX_RAD;
    e->ymax = ymax + MAX_RAD;
    e->zmin = zmin - MAX_RAD;
    e->zmax = zmax + MAX_RAD;
}

/*
 * Radius of a sphere at the origin enclosing every atom.  Assumes the
 * structure has already been centred.
 */
double compute_sphere(const Structure *s)
{
    double rad_max = -999.9;
    int i;

    for (i = 0; i < s->natoms; i++) {
        double x = s->pos[i][X], y = s->pos[i][Y], z = s->pos[i][Z];
        double rad = sqrt(x * x + y * y + z * z) + MAX_RAD;

        if (rad > rad_max)
            rad_max = rad;
    }

    return rad_max;
}

/* ------------------------------------------------------------------ */
/* Bonds                                                              */
/* ------------------------------------------------------------------ */

/*
 * Hydrogen test for the bonding rules below.  This reads the assigned type
 * rather than the atom name, so an element whose symbol merely starts with H
 * -- mercury, hafnium, helium -- is no longer capped at one bond.  Under
 * --legacy-elements the type still comes from the name, so behaviour there is
 * unchanged.
 */
static int is_hydrogen(const Structure *s, int i)
{
    return s->type[i] == H_TYPE;
}

/*
 * Distance-check bond search, with special-casing to stop hydrogens picking
 * up more than one partner.  Originally by M. Pique.
 *
 * The outer index i runs ahead of the inner index j because hydrogens
 * traditionally follow the heavy atom they hang off, so working backwards
 * finds that heavy atom first and the search for i can stop at one bond.
 *
 * Pass bonds == NULL to count without recording; pass a non-NULL bonds to
 * record up to max_bonds pairs.  The 1993 code had two near-copies of this
 * loop -- one to count and one to fill -- which had already drifted apart
 * (one used hypot(), the other an inline sqrt).  A single traversal cannot
 * disagree with itself about how many bonds there are.
 */
int find_bonds(const Structure *s, double threshold, int max_bonds,
               int **bonds)
{
    int i, j;
    int nbonds = 0;

    for (i = s->natoms - 1; i > 0; i--) {
        for (j = i - 1; j >= 0; j--) {
            double dx, dy, dz, dxy;

            /* Never bond hydrogens to each other. */
            if (is_hydrogen(s, i) && is_hydrogen(s, j))
                continue;

            dx = fabs(s->pos[i][X] - s->pos[j][X]);
            if (dx > threshold)
                continue;
            dy = fabs(s->pos[i][Y] - s->pos[j][Y]);
            if (dy > threshold)
                continue;
            dxy = hypot(dx, dy);
            if (dxy > threshold)
                continue;
            dz = fabs(s->pos[i][Z] - s->pos[j][Z]);
            if (dz > threshold)
                continue;
            if (hypot(dxy, dz) > threshold)
                continue;

            if (bonds != NULL) {
                if (nbonds >= max_bonds)
                    return nbonds;
                bonds[nbonds][0] = j;
                bonds[nbonds][1] = i;
            }
            nbonds++;

            if (is_hydrogen(s, i))
                break; /* one bond per hydrogen */
        }
    }

    return nbonds;
}

/* ------------------------------------------------------------------ */
/* POV-Ray output                                                     */
/* ------------------------------------------------------------------ */

/*
 * Derive a legal POV-Ray identifier from the output stem.  POV identifiers
 * are alphanumerics and underscores and may not begin with a digit, so a
 * stem like "out/1crn" -- perfectly good as a path -- cannot be used as one
 * directly.  The old code passed the stem through unfiltered and emitted
 * scenes that would not parse.
 */
static void pov_identifier(const char *stem, char *out, size_t outsz)
{
    const char *base;
    const char *slash;
    size_t i;

    slash = strrchr(stem, '/');
    base = slash ? slash + 1 : stem;

    if (*base == '\0' || isdigit((unsigned char)*base)) {
        /* Leading underscore keeps a numeric PDB code such as 1crn usable. */
        snprintf(out, outsz, "_%.*s", (int)(outsz - 2), base);
    } else {
        snprintf(out, outsz, "%.*s", (int)(outsz - 1), base);
    }

    for (i = 0; out[i] != '\0'; i++) {
        if (!isalnum((unsigned char)out[i]) && out[i] != '_')
            out[i] = '_';
    }
}

/* Texture name for an atom type, or NULL if the type has no texture. */
static const char *atom_object_name(int type, int glass)
{
    switch (type) {
    case H_TYPE:  return glass ? "Atom_Glass_H"  : "Atom_H";
    case N_TYPE:  return glass ? "Atom_Glass_N"  : "Atom_N";
    case C_TYPE:  return glass ? "Atom_Glass_C"  : "Atom_C";
    case P_TYPE:  return glass ? "Atom_Glass_P"  : "Atom_P";
    case O_TYPE:  return glass ? "Atom_Glass_O"  : "Atom_O";
    case S_TYPE:  return glass ? "Atom_Glass_S"  : "Atom_S";
    case CA_TYPE: return glass ? "Atom_Glass_Ca" : "Atom_Ca";
    case FE_TYPE: return glass ? "Atom_Glass_Fe" : "Atom_Fe";
    /*
     * Generic fallback.  Before 2.1 this returned NULL and the atom was
     * dropped without a message, so an ion could disappear from a scene while
     * the header atom count still looked right.
     */
    default:      return glass ? "Atom_Glass_X" : "Atom_X";
    }
}

static const char *radii_include(RadiiSet radii)
{
    switch (radii) {
    case RADII_CPK:      return "atoms_cpk.inc";
    case RADII_COVALENT: return "atoms_covalent.inc";
    case RADII_VDW:
    default:             return "atoms_vdw.inc";
    }
}

/*
 * Camera placement.  The molecule is centred on the origin, so the camera
 * sits back along -z far enough that the widest half-extent subtends
 * CAMERA_HALF_ANGLE_DEG.
 *
 * `right x*image_width/image_height` replaces the old fixed `right <4/3,0,0>`.
 * At 4:3 the two are identical; at any other aspect the fixed form rendered
 * with non-square pixels.  The vertical field of view is unchanged at
 * 2*atan(0.5) = 53.13 degrees, which is the figure quiltwright's docs quote.
 */
static void write_camera(FILE *f, const Extents *e)
{
    double xavg = (e->xmax - e->xmin) / 2.0;
    double yavg = (e->ymax - e->ymin) / 2.0;
    double zavg = (e->zmax - e->zmin) / 2.0;
    double theta = CAMERA_HALF_ANGLE_DEG / DEG_PER_RAD_CAM;
    double dist;

    double distx = -xavg / tan(theta);
    double disty = -yavg / tan(theta);
    double distz = -zavg / tan(theta);

    dist = distx < disty ? distx : disty;
    dist = dist < distz ? dist : distz;

    fprintf(f, "camera {\n");
    fprintf(f, "  location  <0, 0, %.3f>\n", dist);
    fprintf(f, "  look_at   <0, 0, 0>\n");
    fprintf(f, "  direction <0, 0, 1>\n");
    fprintf(f, "  up        <0, 1, 0>\n");
    fprintf(f, "  right     x*image_width/image_height\n");
    fprintf(f, "}\n\n");
}

/*
 * A single light off the upper corner.  In POV-Ray 3.x a light_source is not
 * an object, so the old `object { light_source { ... } }` wrapper is gone.
 */
static void write_light(FILE *f, const Extents *e, int area_light)
{
    double xavg = (e->xmax - e->xmin) / 2.0;
    double theta = LIGHT_HALF_ANGLE_DEG / DEG_PER_RAD_CAM;
    double dist = -xavg / tan(theta);

    fprintf(f, "light_source {\n");
    fprintf(f, "  <%.3f, %.3f, %.3f>\n", e->xmax, e->ymax, dist);
    fprintf(f, "  color White\n");
    if (area_light) {
        fprintf(f, "  area_light <%.3f, 0, 0>, <0, 0, %.3f>, 9, 9\n",
                e->xmax / 2.0, e->xmax / 2.0);
        fprintf(f, "  adaptive 1\n");
        fprintf(f, "  jitter\n");
    }
    fprintf(f, "}\n\n");
}

/*
 * Gradient sky with a cloud shell.  The 2.x colour_map ramps were written as
 * [start, end colourA colourB] pairs; 3.x wants a list of control points, so
 * each ramp becomes its two endpoints.
 */
static void write_sky(FILE *f)
{
    fprintf(f, "// Gradient blue sky with white clouds\n");
    fprintf(f, "sphere {\n");
    fprintf(f, "  <0, 0, 0>, 3000\n");
    fprintf(f, "  pigment {\n");
    fprintf(f, "    gradient y\n");
    fprintf(f, "    color_map {\n");
    fprintf(f, "      [0.0 rgb <0.30, 0.30, 1.00>]\n");
    fprintf(f, "      [0.8 rgb <0.70, 0.70, 1.00>]\n");
    fprintf(f, "      [1.0 rgb <0.90, 0.90, 1.00>]\n");
    fprintf(f, "    }\n");
    fprintf(f, "    scale <3000, 3000, 3000>\n");
    fprintf(f, "    quick_color rgb <0.70, 0.70, 1.00>\n");
    fprintf(f, "  }\n");
    fprintf(f, "  finish {\n");
    fprintf(f, "    ambient 0.7\n");
    fprintf(f, "    diffuse 0  // clouds must not shadow the sky\n");
    fprintf(f, "  }\n");
    fprintf(f, "}\n\n");

    fprintf(f, "sphere {\n");
    fprintf(f, "  <0, 0, 0>, 2590\n");
    fprintf(f, "  pigment {\n");
    fprintf(f, "    bozo\n");
    fprintf(f, "    turbulence 0.5\n");
    fprintf(f, "    color_map {\n");
    fprintf(f, "      [0.0 rgbf <1.00, 1.00, 1.00, 1.00>]\n");
    fprintf(f, "      [0.6 rgbf <1.00, 1.00, 1.00, 1.00>]\n");
    fprintf(f, "      [0.8 rgb  <1.00, 1.00, 1.00>]\n");
    fprintf(f, "      [1.0 rgb  <0.80, 0.80, 0.80>]\n");
    fprintf(f, "    }\n");
    fprintf(f, "    quick_color rgb <0.70, 0.70, 1.00>\n");
    fprintf(f, "    scale <1000, 200, 1000>\n");
    fprintf(f, "  }\n");
    fprintf(f, "  finish { ambient 0.7 diffuse 0 }\n");
    fprintf(f, "}\n\n");
}

static void write_sky_plain(FILE *f)
{
    fprintf(f, "// Plain white sky\n");
    fprintf(f, "sphere {\n");
    fprintf(f, "  <0, 0, 0>, 3000\n");
    fprintf(f, "  pigment { color White }\n");
    fprintf(f, "  finish { ambient 0.7 diffuse 0 }\n");
    fprintf(f, "}\n\n");
}

/* Y of the ground plane, set below the lowest atom. */
static double ground_level(const Extents *e)
{
    double yavg = (e->ymax - e->ymin) / 2.0;

    return -1.0 * (fabs(e->ymin) + fabs(yavg) / 2.0);
}

static void write_ground(FILE *f, const Extents *e)
{
    fprintf(f, "plane {\n");
    fprintf(f, "  y, %.3f\n", ground_level(e));
    fprintf(f, "  pigment { color RichBlue quick_color RichBlue }\n");
    fprintf(f, "  normal { bumps 0.2 }\n");
    fprintf(f, "}\n\n");
}

static void write_check(FILE *f, const Extents *e)
{
    double xavg = (e->xmax - e->xmin) / 2.0;
    double yavg = (e->ymax - e->ymin) / 2.0;

    fprintf(f, "// The beloved raytracer checkerboard\n");
    fprintf(f, "plane {\n");
    fprintf(f, "  y, %.3f\n", ground_level(e));
    fprintf(f, "  pigment {\n");
    fprintf(f, "    checker color Black color White\n");
    fprintf(f, "    scale %.3f\n", (xavg + yavg) / 3.0);
    fprintf(f, "  }\n");
    fprintf(f, "  finish { ambient 0.2 diffuse 0.8 }\n");
    fprintf(f, "}\n\n");
}

/*
 * Emit one merge/union member per atom.
 *
 * Under --legacy-elements an atom with no dedicated texture is dropped, which
 * is what happened before 2.1; reproducing an old render means reproducing
 * the omission too.  It is counted and reported either way, so the drop is no
 * longer silent.  By default such atoms render as the generic Atom_X.
 */
static void write_atoms(FILE *f, const Structure *s, const Options *opt,
                        const char *scale_id, int glass)
{
    int i;

    for (i = 0; i < s->natoms; i++) {
        const char *name;

        if (opt->legacy_elements && s->type[i] == ANY_TYPE)
            continue;

        name = atom_object_name(s->type[i], glass);
        if (name == NULL)
            continue;

        fprintf(f,
                "  object { %s scale <%s, %s, %s> translate <%.3f, %.3f, "
                "%.3f> }\n",
                name, scale_id, scale_id, scale_id, s->pos[i][X],
                s->pos[i][Y], s->pos[i][Z]);
    }
}

static void write_header(FILE *f, const Structure *s, const Extents *e,
                         double sphere_rad)
{
    char stamp[32];
    time_t now = time(NULL);
    struct tm *lt = localtime(&now);

    /*
     * strftime replaces a hand-rolled sprintf that wrote "%02d/%02d/%02d" from
     * tm_year -- years since 1900 -- into a 10-byte buffer.  In 2026 that
     * printed the year as "126" and overran the buffer, which is why the
     * build had to disable _FORTIFY_SOURCE.
     */
    if (lt == NULL || strftime(stamp, sizeof stamp, "%Y-%m-%d %H:%M:%S", lt) == 0)
        snprintf(stamp, sizeof stamp, "unknown");

    fprintf(f, "//\n");
    fprintf(f, "// Prepared by pdb2pov %s on %s\n", PDB2POV_VERSION, stamp);
    fprintf(f, "// Author: Eric G. Suchanek, Ph.D.\n");
    fprintf(f, "//\n");
    fprintf(f, "//\tAtoms: %4d\n", s->natoms);
    fprintf(f, "//\tExtent:\tXmin: %.3f Xmax: %.3f,\n//\t\tYmin: %.3f, Ymax: %.3f\n",
            e->xmin, e->xmax, e->ymin, e->ymax);
    fprintf(f, "//\t\tZmin: %.3f Zmax: %.3f\n", e->zmin, e->zmax);
    fprintf(f, "//\tEnclosing Sphere: %.3f\n", sphere_rad);
    fprintf(f, "//\n\n");
}

int write_output(const Options *opt, const Structure *s, int **bonds,
                 int nbonds)
{
    char path[PATH_MAX_LEN];
    char ident[PATH_MAX_LEN];
    Extents e;
    FILE *f;
    double sphere_rad;
    double atom_scale = opt->radii_scale;
    int i;

    compute_extents(s, &e);
    sphere_rad = compute_sphere(s);
    sphere_rad += SPHERE_FUDGE * sphere_rad;

    /* Ball-and-stick shrinks the spheres so the bond cylinders show. */
    if (opt->ball_stick)
        atom_scale *= 0.3;

    if (!join_extension(path, sizeof path, opt->output_stem,
                        opt->object_only ? ".inc" : ".pov"))
        return ERR_CANT_WRITE_OUTPUT;

    pov_identifier(opt->output_stem, ident, sizeof ident);

    f = fopen(path, "w");
    if (f == NULL) {
        fprintf(stderr, "pdb2pov: can't open '%s' for writing\n", path);
        return ERR_CANT_WRITE_OUTPUT;
    }

    printf("Writing output file <%s>\n", path);

    write_header(f, s, &e, sphere_rad);

    if (opt->object_only) {
        /*
         * An include has to leave the language version as it found it, so the
         * host scene is not silently switched to 3.7 from whatever it set.
         */
        fprintf(f, "#declare %s_pov_version = version;\n", ident);
        fprintf(f, "#version %s;\n\n", POV_VERSION);
    } else {
        fprintf(f, "#version %s;\n\n", POV_VERSION);
        /* Required from 3.7 onwards; POV-Ray errors out without it. */
        fprintf(f, "global_settings { assumed_gamma 1.0 }\n\n");

        fprintf(f, "#include \"colors.inc\"\n");
        fprintf(f, "#include \"shapes.inc\"\n");
        fprintf(f, "#include \"textures.inc\"\n");
        fprintf(f, "#include \"%s\"\n", radii_include(opt->radii));
        fprintf(f, "#include \"atoms2.inc\"\n");
        if (opt->glass_atoms)
            fprintf(f, "#include \"atoms_glass2.inc\"\n");
        fprintf(f, "\n");

        write_camera(f, &e);
        write_light(f, &e, opt->area_light);

        if (opt->backdrop == BACKDROP_SKY)
            write_sky(f);
        else if (opt->backdrop == BACKDROP_PLAIN)
            write_sky_plain(f);

        if (opt->ground == GROUND_PLAIN)
            write_ground(f, &e);
        else if (opt->ground == GROUND_CHECKER)
            write_check(f, &e);
    }

    /*
     * Float declarations need a terminating semicolon from POV-Ray 3.5
     * onwards.  Object and texture declarations do not, and do not get one.
     */
    if (opt->ball_stick)
        fprintf(f, "#declare BOND_RADIUS = %.2f;\n", BOND_RAD);
    if (opt->glass_atoms)
        fprintf(f, "#declare ATM_SCL_B = %.2f;\n", atom_scale / 0.3);
    fprintf(f, "#declare ATM_SCL = %.2f;\n\n", atom_scale);

    /* The glass shell, if asked for, is merged so interior faces vanish. */
    if (opt->glass_atoms) {
        fprintf(f, "#declare %s_obj_glass = merge {\n", ident);
        write_atoms(f, s, opt, "ATM_SCL_B", 1);
        fprintf(f, "}\n\n");
    }

    fprintf(f, "#declare %s_obj = union {\n", ident);
    write_atoms(f, s, opt, "ATM_SCL", 0);

    if (opt->ball_stick && bonds != NULL) {
        for (i = 0; i < nbonds; i++) {
            int a = bonds[i][0];
            int b = bonds[i][1];

            fprintf(f,
                    "  cylinder { <%.3f, %.3f, %.3f>, <%.3f, %.3f, %.3f>, "
                    "BOND_RADIUS texture { White_Bond } }\n",
                    s->pos[a][X], s->pos[a][Y], s->pos[a][Z], s->pos[b][X],
                    s->pos[b][Y], s->pos[b][Z]);
        }
    }

    if (opt->glass_atoms)
        fprintf(f, "  object { %s_obj_glass }\n", ident);

    fprintf(f, "}\n\n");

    /*
     * The enclosing sphere radius, as a float a host scene can actually read.
     * The 1993 code spent it on `bounded_by { sphere { ... } }`, which
     * POV-Ray 3.x rejects as redundant -- it warns "CSG union unnecessarily
     * bounded", because its automatic bounding is both tighter than a sphere
     * around the whole molecule and guaranteed not to clip.  The number is
     * still worth publishing: it is what sets the depth budget when these
     * scenes are rendered for a holographic display.
     */
    fprintf(f, "#declare %s_enclosing_radius = %.3f;\n", ident, sphere_rad);
    fprintf(f, "#declare %s = object { %s_obj }\n", ident, ident);

    if (opt->object_only)
        fprintf(f, "\n#version %s_pov_version;\n", ident);
    else
        fprintf(f, "\nobject { %s }\n", ident);

    if (fclose(f) != 0) {
        fprintf(stderr, "pdb2pov: error writing '%s'\n", path);
        return ERR_CANT_WRITE_OUTPUT;
    }

    return 0;
}

/* ------------------------------------------------------------------ */
/* Structure lifetime                                                 */
/* ------------------------------------------------------------------ */

static void structure_free(Structure *s, int capacity)
{
    free_ivector(s->type);
    free_cmatrix(s->element, (size_t)capacity);
    free_cmatrix(s->res_name, (size_t)capacity);
    free_cmatrix(s->atom_name, (size_t)capacity);
    free_dvector(s->charge);
    free_dmatrix(s->pos, (size_t)capacity);
    memset(s, 0, sizeof *s);
}

static int structure_alloc(Structure *s, int capacity)
{
    memset(s, 0, sizeof *s);

    s->pos = dmatrix((size_t)capacity, 3);
    s->charge = dvector((size_t)capacity);
    s->atom_name = cmatrix((size_t)capacity, ATOM_NAME_LEN);
    s->res_name = cmatrix((size_t)capacity, RES_NAME_LEN);
    s->element = cmatrix((size_t)capacity, ELEM_NAME_LEN);
    s->type = ivector((size_t)capacity);

    if (s->pos == NULL || s->charge == NULL || s->atom_name == NULL ||
        s->res_name == NULL || s->element == NULL || s->type == NULL) {
        structure_free(s, capacity);
        return 0;
    }

    return 1;
}

/* ------------------------------------------------------------------ */
/* Entry point                                                        */
/* ------------------------------------------------------------------ */

int main(int argc, char **argv)
{
    Options opt;
    Structure s;
    ParseStats stats;
    double center[3];
    int **bonds = NULL;
    int nbonds = 0;
    int capacity;
    int status;

    if (!parse_options(argc, argv, &opt))
        return -ERR_PARSE_ARGS;

    memset(&stats, 0, sizeof stats);

    printf("Scanning atom file <%s>... ", opt.input_path);
    fflush(stdout);

    capacity = count_pdb_atoms(opt.input_path);
    if (capacity < 0)
        return -ERR_CANT_READ_INPUT;
    if (capacity == 0) {
        printf("\npdb2pov: no atoms found in <%s>\n", opt.input_path);
        return -ERR_NO_ATOMS;
    }

    printf("got <%d> atoms.\n", capacity);

    if (!structure_alloc(&s, capacity)) {
        fprintf(stderr, "pdb2pov: out of memory for %d atoms\n", capacity);
        return -ERR_NO_ATOMS;
    }

    if (opt.atm_format)
        status = read_atm(opt.input_path, capacity, &s, &stats);
    else
        status = read_pdb(opt.input_path, capacity, &opt, &s, &stats);

    if (status < 0 || s.natoms == 0) {
        fflush(stdout); /* keep the diagnostic after the progress lines */
        if (stats.skipped_chain > 0 && stats.accepted == 0)
            fprintf(stderr,
                    "pdb2pov: no atoms left after filtering to chain(s) "
                    "'%s'\n", opt.chain_filter);
        else
            fprintf(stderr, "pdb2pov: couldn't read atoms from <%s>\n",
                    opt.input_path);
        structure_free(&s, capacity);
        return -ERR_NO_ATOMS;
    }

    make_atom_types(&opt, &s, &stats);

    if (s.natoms != capacity)
        printf("Kept <%d> of <%d> atoms.\n", s.natoms, capacity);
    report_parse(&stats, opt.legacy_elements);

    printf("Re-orienting and positioning structure.\n");
    rotate_structure(&s, opt.xrot, opt.yrot, opt.zrot);
    center_structure(&s, center);
    flip_zaxis(&s);

    if (opt.ball_stick) {
        printf("Computing bonds: ");
        fflush(stdout);

        nbonds = find_bonds(&s, opt.bond_threshold, 0, NULL);
        printf("found %d bonds.\n", nbonds);

        if (nbonds <= 0) {
            fprintf(stderr,
                    "pdb2pov: no bonds within %.2f angstroms; try a larger "
                    "-d threshold\n",
                    opt.bond_threshold);
            structure_free(&s, capacity);
            return -ERR_NO_BONDS;
        }

        bonds = imatrix((size_t)nbonds, 2);
        if (bonds == NULL) {
            fprintf(stderr, "pdb2pov: out of memory for %d bonds\n", nbonds);
            structure_free(&s, capacity);
            return -ERR_CANT_ALLOC_BONDS;
        }

        find_bonds(&s, opt.bond_threshold, nbonds, bonds);
    }

    status = write_output(&opt, &s, bonds, nbonds);

    free_imatrix(bonds, (size_t)nbonds);
    structure_free(&s, capacity);

    return status == 0 ? 0 : -status;
}

/* ------------------------- End of pdb2pov.c ------------------------- */
