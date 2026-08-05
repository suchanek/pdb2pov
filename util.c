/*
 * util.c -- rectangular array allocators for pdb2pov.
 *
 * Copyright (c) 1993-2026 Eric G. Suchanek, Ph.D.  All rights reserved.
 * Subject to the GNU License.
 *
 * History: these began as the Numerical-Recipes dmatrix/dvector family.  Three
 * things about that lineage are gone:
 *
 *   - the returned pointer was offset by a subscript base (`- nrl`) that was
 *     always zero.  Computing it was undefined behaviour for any non-zero
 *     base, and it bought nothing at zero.
 *   - the row-pointer array was one entry longer than the row count, because
 *     the bounds were inclusive high subscripts rather than counts.
 *   - an allocation failure called exit() from inside the allocator, so a
 *     caller could not release what it already held.  They now return NULL.
 */

#include <stdlib.h>

#include "util.h"

/*
 * Allocating a rows x cols array is the same shape of work for each element
 * type, so the three matrix allocators share one implementation.  The row
 * pointers are stored in a void * array, which is layout-compatible with the
 * T * array the caller sees.
 */
static void **alloc_matrix(size_t rows, size_t cols, size_t elem_size)
{
    void **m;
    size_t i;

    if (rows == 0 || cols == 0)
        return NULL;

    /* Guard the row-byte count against overflow before asking for it. */
    if (cols > (size_t)-1 / elem_size)
        return NULL;

    m = calloc(rows, sizeof *m);
    if (m == NULL)
        return NULL;

    for (i = 0; i < rows; i++) {
        m[i] = calloc(cols, elem_size);
        if (m[i] == NULL) {
            /* calloc() zeroed the tail, so freeing the whole array is safe. */
            while (i > 0)
                free(m[--i]);
            free(m);
            return NULL;
        }
    }

    return m;
}

static void free_matrix(void **m, size_t rows)
{
    size_t i;

    if (m == NULL)
        return;

    for (i = 0; i < rows; i++)
        free(m[i]);
    free(m);
}

double **dmatrix(size_t rows, size_t cols)
{
    return (double **)alloc_matrix(rows, cols, sizeof(double));
}

void free_dmatrix(double **m, size_t rows)
{
    free_matrix((void **)m, rows);
}

int **imatrix(size_t rows, size_t cols)
{
    return (int **)alloc_matrix(rows, cols, sizeof(int));
}

void free_imatrix(int **m, size_t rows)
{
    free_matrix((void **)m, rows);
}

char **cmatrix(size_t rows, size_t cols)
{
    return (char **)alloc_matrix(rows, cols, sizeof(char));
}

void free_cmatrix(char **m, size_t rows)
{
    free_matrix((void **)m, rows);
}

double *dvector(size_t n)
{
    return n ? calloc(n, sizeof(double)) : NULL;
}

void free_dvector(double *v)
{
    free(v);
}

int *ivector(size_t n)
{
    return n ? calloc(n, sizeof(int)) : NULL;
}

void free_ivector(int *v)
{
    free(v);
}
