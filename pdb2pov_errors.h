/* $Id: pdb2pov_errors.h,v 1.1 1993/11/25 19:59:54 eric Exp eric $ */

/* $Log: pdb2pov_errors.h,v $
 * Revision 1.1  1993/11/25  19:59:54  eric
 * Initial revision
 * */


#define ERR_PARSE_ARGS -2
#define ERR_NO_ATOMS -3
#define ERR_NO_BONDS -4
#define ERR_CANT_WRITE_OUTPUT -5
#define ERR_CANT_READ_INPUT -6
#define ERR_READ_INPUT -7
#define ERR_CANT_ALLOC_BONDS -8

static char *pdb2pov_errorstringsp[] = 
{0, 0, "Can't parse arguments", "Couldn't read any atoms",
	"Couldn't compute any bonds", "Couldn't write output file",
	"Couldn't read input file", "Couldn't allocate memory for bond array"};
