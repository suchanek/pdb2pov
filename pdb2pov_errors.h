/*
 * pdb2pov_errors.h -- failure codes.
 *
 * The values are negative for historical reasons.  main() negates them on the
 * way out, because a process exit status is an unsigned byte and returning a
 * negative value from main() reports as its two's complement.
 *
 * Copyright (c) 1993-2026 Eric G. Suchanek, Ph.D.  All rights reserved.
 * Subject to the GNU License.
 */

#ifndef PDB2POV_ERRORS_H
#define PDB2POV_ERRORS_H

#define ERR_PARSE_ARGS        -2
#define ERR_NO_ATOMS          -3
#define ERR_NO_BONDS          -4
#define ERR_CANT_WRITE_OUTPUT -5
#define ERR_CANT_READ_INPUT   -6
#define ERR_READ_INPUT        -7
#define ERR_CANT_ALLOC_BONDS  -8

#endif /* PDB2POV_ERRORS_H */
