/* Copyright (c) 1992 SAS Institute, Inc., Cary, NC USA */
/* All Rights Reserved */

#ifdef _M68881
double __builtin_fpc(int, double);

/*  for optimal 040 use of this file, modify math.h as follows

#ifdef _M68040
#define __FUNC_CHARS __inline
#else
#define __FUNC_CHARS
#endif

and for each of the following functions put a __FUNC_CHARS in the
definiton in math.h, for example:

double __FUNC_CHARS sin(double);

put a __FUNC_CHARS in front of the following functions in math.h
    pow
    floor
    fmod
    acos
    asin
    atan
    atan2
    cos
    exp
    log
    log10
    sin
    tan
*/

#ifdef _M68040
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#define ___ERRNUM HUGE_VAL
#endif

#ifndef __FUNC_CHARS /* this would mean math.h has not been modified */
#define __FUNC_CHARS 
#endif

/* need to move sqrt defn to top since it is used in funcs for 040 */
#define sqrt(d)     __builtin_fpc(0x0004,d)

#ifndef _M68040
#define acos(d)     __builtin_fpc(0x001C,d)
#else
#define __m68pid2 1.5707963267948966192

double __FUNC_CHARS acos(x)
double x;
{
     return(__m68pid2-asin(x));
}
#endif
#ifndef _M68040
#define asin(d)     __builtin_fpc(0x000C,d)
#else

double __FUNC_CHARS asin(x)
double x;
{
     double y,y2,sum;
     int sflag, yflag;
     if ( x < 0 ) { sflag = 1; x = -x; }
     else sflag = 0;
     if ( x > 0.5 ) {
          if ( x > 1.0 ) {
               if ( sflag ) x = -x;
               return(___ERRNUM);
          }
          y = sqrt(0.5*(1.-x));
          yflag = 1;
     }
     else { y = x; yflag = 0; }
     y2 = y*y;
     sum = y * ( 1.00000000000000001365+y2*(
                   0.16666666666664838839+y2*(
                   0.07500000000404881646+y2*(
                   0.04464285679119552341+y2*(
                   0.03038196027237537659+y2*(
                   0.02237173658111094020+y2*(
                   0.01735997031086759961+y2*(
                   0.01388302431756389624+y2*(
                   0.01218170270685370993+y2*(
                   0.00648075111416853345+y2*(
                   0.01964268701248804731+y2*(
                  -0.01638489913725409907+y2*(
                   0.03201166198466493988)))))))))))));
     if ( yflag ) sum = __m68pid2 - 2.*sum;
     if ( sflag ) return(-sum);
     else return(sum);
}
#endif
#ifndef _M68040
#define atan(d)     __builtin_fpc(0x000A,d)
#else
double __FUNC_CHARS atan2(y,x)
double y,x;
{
/*
     arctangent between -pi and pi.
*/
     int sflag, piflag, pflag, qflag;
     double z,sum,z2;
     if ( y < 0 ) { y = -y;  sflag = 1; }
     else  sflag = 0;
     if ( x < 0 ) { x = -x; piflag = 1; }
     else piflag = 0;
     if ( y > x ) { pflag = 1; z = y; y = x; x = z; }
     else pflag = 0;
/*
     We are now in the range 0-1
*/
     if ( y == 0 && x == 0 ) {
          return(0);
     }
     if ( y > 0.2679491924311227065*x ) {
          z = (y*1.7320508075688772935-x)/(y+1.7320508075688772935*x);
          qflag = 1;
     }
     else {
          z = y/x;
          qflag = 0;
     }
/*
     Now we are in the range -(2-sqrt(3))<z<2-sqrt(3)
*/
     z2 = z*z;
     sum = z*( 0.99999999999999998493+z2*(
             -0.33333333333329920123+z2*(
              0.19999999998726326765+z2*(
             -0.14285714102479574031+z2*(
              0.11111097879199028445+z2*(
             -0.09090370583742286501+z2*(
              0.07679359531984028197+z2*(
             -0.06483115238224471189+z2*(
              0.04438679771995954480)))))))));
     if ( qflag ) sum += 0.5235987755982988731;
     if ( pflag ) sum = 1.5707963267948966192 - sum;
     if ( piflag ) sum = 3.1415926535897932384 - sum;
     if ( sflag ) sum = -sum;
     return(sum);
}

double __FUNC_CHARS atan(x)
double x;
{
     return(atan2(x,1.));
}
#endif
#ifndef _M68040
#define cos(d)      __builtin_fpc(0x001D,d)
#else
#define __m68pid4 0.7853981633974483096
#define __m68pi   3.1415926535897932384
#define __m68pi2  6.2831853071795864769

double __FUNC_CHARS floor(x)
double x;
{
/*
     Routine determines the largest integer that is smaller than x.
     The return value is a double again.
     We have to use the IEEE explicitly here.
     Routine made by J.A.M. Vermaseren 20-apr-1992
*/
     double scr, y;
     int n;
     y = x;
     if ( y < 0 ) y = -y;
     scr = y;
     n = *((short *)(&scr));
     if ( n < 0 ) n = -n;
     n -= 0x3FF0;
     n >>= 4;
     if ( n < 0 ) { /* exponent < 0 --> trivial */
          if ( x < 0 ) return(-1.);
          else return(0.);
     }
/*
     There are 53 bits of which 52 after the period. This means that
     if n >= 52 we cannot do anything.
*/
     if ( n >= 52 ) return(x);
     if ( n <= 20 ) {
          ((long *)(&scr))[1] = 0;
          n = 20 - n;
          ((long *)(&scr))[0] &= (-1) << n;
     }
     else {
          n = 52 - n;
          ((long *)(&scr))[1] &= (-1) << n;
     }
     if ( x < 0 ) x = -scr-1;
     else x = scr;
     return(x);
}

double __FUNC_CHARS fmod(x,y)
double x,y;
{
/*
     Routine determines the remainder of x/y with the same sign as x
     Routine made by J.A.M. Vermaseren 20-apr-1992
*/
     double z;
     if ( y == 0. ) return(x);
     if ( y > 0 ) {
          if ( x < y && x > -y ) return(x);
     }
     else if ( x > y && x < -y ) return(x);
     z = x/y;
     z -= floor(z);
     z *= y;
     if ( x < 0 ) {
          if ( z > 0 ) {
               if ( y > 0 ) z = z - y;
               else z = z + y;
          }
     }
     else {
          if ( z < 0 ) {
               if ( y > 0 ) z = z + y;
               else z = z - y;
          }
     }
     return(z);
}


double __FUNC_CHARS cos(x)
double x;
{
     double x2,sum;
     int sflag;
     if ( x < 0 ) x = -x;
     if ( x >= __m68pi2 ) x = fmod(x,__m68pi2);
     if ( x > __m68pi ) x = __m68pi2-x;
     if ( x > __m68pid2 ) { sflag = 1; x = __m68pi - x; }
     else sflag = 0;
/*
     Now we are in the first quadrant
*/
     if ( x < __m68pid4 ) {
          x2 = x*x;
          sum =  0.99999999999999995287+x2*(
              -0.49999999999999251135+x2*(
               0.04166666666647237368+x2*(
              -0.00138888888699833768+x2*(
               0.00002480157854004439+x2*(
              -0.00000027555234099211+x2*(
               0.00000000206304660116))))));
     }
     else {
          x = __m68pid2 - x;
          x2 = x*x;
          sum = x*( 0.99999999999999999685+x2*(
                 -0.16666666666666616691+x2*(
                  0.00833333333332036780+x2*(
                 -0.00019841269828654287+x2*(
                  0.00000275573133777316+x2*(
                 -0.00000002505071717794+x2*(
                  0.00000000015894748309)))))));
     }

     if ( sflag ) return(-sum);
     else return(sum);
}
#endif
#define cosh(d)     __builtin_fpc(0x0019,d)
#ifndef _M68040
#define exp(d)      __builtin_fpc(0x0010,d)
#else
#define __m68log2  0.6931471805599453094
#define __m68ilog2  (1./0.6931471805599453094)

double __FUNC_CHARS exp(x)
double x;
{
/*
     routine for the evaluation of the exp (IEEE)
     Method:
          1: y = x/log(2)  -> exp(x) = 2^y;
          2: z = n+y, -0.5 <= y < 0.5
          3: exp(x) = exp(z*log(2))*2^n
     For the last series we need only terms till 13!
     Routine made by J.A.M.Vermaseren 20-apr-1992
*/
     double y,z,scr;
     long n;
     int m,sflag;
     y = x*__m68ilog2;
     if ( y < 0 ) { sflag = 1; y = -y; }
     else sflag = 0;
     scr = y + 0.5;
     n = *((long *)(&scr));
     m = (n >> 20) - 0x3FF;
     if ( m > 7 ) {
          return(___ERRNUM);
     }
     n &= 0xFFFFF;
     n |= 0x100000;
     n >>= (20-m);
     z = (y - n)*__m68log2;
     if ( sflag ) { n = -n; z = -z; }
     scr = 0.99999999999999999693+z*(
          0.99999999999999999850+z*(
          0.50000000000000183922+z*(
          0.16666666666666701804+z*(
          0.04166666666648804594+z*(
          0.00833333333330984930+z*(
          0.00138888889523274996+z*(
          0.00019841269908532560+z*(
          0.00002480148546912669+z*(
          0.00000275572255877329+z*(
          0.00000027632644397330+z*(
          0.00000002511466084085)))))))))));
     n <<= 20;
     *((long *)(&scr)) += n;
     return(scr);
}
#endif
#define fabs(d)     __builtin_fpc(0x0018,d)
#ifndef _M68040
/* #define floor(d)    __builtin_fpc(0x0003,d) */ /* works only on pos numbers */
#endif

#ifndef _M68040
#define log(d)      __builtin_fpc(0x0014,d)
#else
double __FUNC_CHARS log(x)
double x;
{
/*
     Routine for log(x) in IEEE double precision.
     Note that we have to hack around in the exponent!
     Steps: 1 x = y*2^n, with y between sqrt(2) and 1/sqrt(2)
            2 convert to z = (y-1)/(y+1)
            3 ln = 2*z*sum(i=0,9, z^(2*i)/(2*i+1)) + n*ln(2)
     Routine made by J.A.M. Vermaseren 20-apr-1992.
*/
     double scr, y, z, z2, sum;
     short n;
     if ( x < 0 ) {
          return(___ERRNUM);
     }
     scr = x;
     n = ( (*((short *)(&scr))) & 0x7FF0 ) - 0x3FF0;
     *((short *)(&scr)) -= n;
     n >>= 4;
     if ( scr > 1.4142135623730950488 ) {
          *((short *)(&scr)) -= 0x10;
          n++;
     }
     y = scr;
     z = (y-1.)/(y+1.);
     z2 = z*z;
     sum = 2.*z*(0.99999999999999999886+z2*(
                0.33333333333333825805+z2*(
                0.19999999999649522967+z2*(
                0.14285714380662481221+z2*(
                0.11111098494595009814+z2*(
                0.09091817553294820737+z2*(
                0.07656220982777783836+z2*(
                0.07405280893003482131))))))))
           +n*0.6931471805599453094;
     return(sum);
}
#endif
#ifndef _M68040
#define log10(d     __builtin_fpc(0x0015,d)
#else
double __FUNC_CHARS log10(x)
double x;
{
     return(log(x)*0.43429448190325182765);
}
#endif
#ifndef _M68040
#define sin(d)      __builtin_fpc(0x000E,d)
#else

double __FUNC_CHARS sin(x)
double x;
{
     double x2,sum;
     int sflag;
     if ( x < 0 ) { sflag = 1; x = -x; }
     else sflag = 0;
     if ( x >= __m68pi2 ) x = fmod(x,__m68pi2);
     if ( x >= __m68pi ) { x = __m68pi2 - x; sflag = 1-sflag; }
     if ( x > __m68pid2 ) x = __m68pi - x;
/*
     Now we are in the first quadrant
*/
     if ( x > __m68pid4 ) {
          x = __m68pid2 - x;
          x2 = x*x;
          sum =  0.99999999999999995287+x2*(
              -0.49999999999999251135+x2*(
               0.04166666666647237368+x2*(
              -0.00138888888699833768+x2*(
               0.00002480157854004439+x2*(
              -0.00000027555234099211+x2*(
               0.00000000206304660116))))));
     }
     else {
          x2 = x*x;
          sum = x*( 0.99999999999999999685+x2*(
                 -0.16666666666666616691+x2*(
                  0.00833333333332036780+x2*(
                 -0.00019841269828654287+x2*(
                  0.00000275573133777316+x2*(
                 -0.00000002505071717794+x2*(
                  0.00000000015894748309)))))));
     }

     if ( sflag ) return(-sum);
     else return(sum);
}
#endif
#define sinh(d)     __builtin_fpc(0x0002,d)
#ifndef _M68040
#define tan(d)      __builtin_fpc(0x000F,d)
#else

double __FUNC_CHARS tan(x)
double x;
{
/*
     Routine for the tangent. It puts x in the range 0-pi/4.
     Then it uses a 6,6 Maehly approximation.
     Derived and coded by J.A.M.Vermaseren 21-apr-1992.
*/

     double sum,x2;
     int sflag, pflag;

     if ( x < 0 ) { sflag = 1; x = -x; }
     else sflag = 0;
     if ( x >= __m68pi ) x = fmod(x,__m68pi);
     if ( x >= __m68pid2 ) { x = __m68pi - x; sflag = 1 - sflag; }
     if ( x > __m68pid4 ) { pflag = 1; x = __m68pid2 - x; }
     else pflag = 0;
     x2 = x*x;
     sum = x*(-4797.60570220792587597407+x2*(
               615.45028715821815261900+x2*(
               -13.46133040142504243340+x2*
                 0.03590005535592163977)))/(
             -4797.60570220792579574817+x2*(
              2214.65218789418096945538+x2*(
              -111.99796607144895992639+x2)));
     if ( pflag ) sum = 1./sum;
     if ( sflag ) return(-sum);
     else return(sum);
}
#endif
#define tanh(d)     __builtin_fpc(0x0009,d)

#define fintrz(d)   __builtin_fpc(0x0003,d)
#ifndef _M68040
#define flognp1(d)  __builtin_fpc(0x0006,d)
#else
#define flognp1(d)  log((d)+1.0)
#endif
#ifndef _M68040
#define fetoxm1(d)  __builtin_fpc(0x0008,d)
#else
#define fetoxm1(d)  (exp((d)) - 1.0)
#endif
#define fatanh(d)   __builtin_fpc(0x000D,d)
#ifndef _M68040
#define fetox(d)    __builtin_fpc(0x0010,d)
#else
#define fetox exp
#endif
#ifdef _M68040
double __FUNC_CHARS pow(x,y)
double x,y;
{
/*
     Forget problems with pow(0,0)
*/
     if ( x == 0. ) return(0.);
     if ( y == 0. ) return(1.);
     return(exp(y*log(x)));
}
#endif
#ifndef _M68040
#define ftwotox(d)  __builtin_fpc(0x0011,d)
#else
#define ftwotox(d) pow(2.0,d)
#endif
#ifndef _M68040
#define ftentox(d)  __builtin_fpc(0x0012,d)
#else
#define ftentox(d) pow(10.0,d)
#endif
#ifndef _M68040
#define flogn(d)    __builtin_fpc(0x0014,d)
#else
#define flogn log
#endif
#define fneg(d)     __builtin_fpc(0x001A,d)
#define fgetman(d)  __builtin_fpc(0x001F,d)
#endif
