/*
 * pdb2pov_protos.h -- prototypes for the allocation routines in util.c.
 *
 * pdb2pov.c calls dmatrix(), dvector(), cmatrix(), ivector() and imatrix()
 * without declaring them.  Under K&R rules an undeclared function is assumed
 * to return int, so on a 64-bit host the pointers they return are truncated
 * to 32 bits and the program segfaults on first use.  This was harmless when
 * the code was written and is fatal now.
 *
 * Force-included from the Makefile (see PORTFLAGS) so the 1993 sources need
 * no edits.
 *
 * Copyright (c) 1993-2026 Eric G. Suchanek, Ph.D.  All rights reserved.
 */

#ifndef PDB2POV_PROTOS_H
#define PDB2POV_PROTOS_H

double **dmatrix(int nrh, int nch);
void     free_dmatrix(double **m, int nrh, int nch);

int    **imatrix(int nrh, int nch);
void     free_imatrix(int **m, int nrh, int nch);

char   **cmatrix(int nrh, int nch);
void     free_cmatrix(char **m, int nrh, int nch);

double  *dvector(int nh);
void     free_dvector(double *v, int nh);

int     *ivector(int nh);
void     free_ivector(int *v, int nh);

#endif /* PDB2POV_PROTOS_H */
