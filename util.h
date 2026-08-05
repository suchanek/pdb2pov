/*
 * util.h -- rectangular array allocators for pdb2pov.
 *
 * These replace the Numerical-Recipes-style routines the 1993 code used.  The
 * originals took (nrh, nch) high-subscript bounds, allocated one extra row,
 * offset every returned pointer by a zero-valued nrl, and called exit() from
 * inside the allocator on failure.  These take plain row/column counts and
 * return NULL so the caller can clean up.
 *
 * Copyright (c) 1993-2026 Eric G. Suchanek, Ph.D.  All rights reserved.
 * Subject to the GNU License.
 */

#ifndef UTIL_H
#define UTIL_H

#include <stddef.h>

/*
 * Each allocator returns a rows x cols array, or NULL if any allocation fails
 * (in which case nothing is left allocated).  The matching free() accepts
 * NULL and tolerates being called on a partially built array.
 */

double **dmatrix(size_t rows, size_t cols);
void free_dmatrix(double **m, size_t rows);

int **imatrix(size_t rows, size_t cols);
void free_imatrix(int **m, size_t rows);

char **cmatrix(size_t rows, size_t cols);
void free_cmatrix(char **m, size_t rows);

double *dvector(size_t n);
void free_dvector(double *v);

int *ivector(size_t n);
void free_ivector(int *v);

#endif /* UTIL_H */
