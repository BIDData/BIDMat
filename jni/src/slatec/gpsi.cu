/* psi.f -- translated by f2c (version 20100827).
   You must link the resulting object file with libf2c:
	on Microsoft Windows system, link with libf2c.lib;
	on Linux or Unix systems, link with .../path/to/libf2c.a -lm
	or, if you install libf2c.a in a standard place, with -lf2c -lm
	-- in that order, at the end of the command line, as in
		cc *.o -lf2c -lm
	Source for libf2c is in /netlib/f2c/libf2c.zip, e.g.,

		http://www.netlib.org/f2c/libf2c.zip
*/

#include "f2c.h"

/* Table of constant values */

extern __device__ integer c__3;
extern __device__ integer c__23;
extern __device__ integer c__16;
extern __device__ integer c__4;
extern __device__ integer c__2;
extern __device__ integer c__1;

__device__ __constant__ real psics[23] = { -.038057080835217922f,.49141539302938713f,
	    -.05681574782124473f,.008357821225914313f,-.001333232857994342f,
	    2.20313287069308e-4f,-3.7040238178456e-5f,6.283793654854e-6f,
	    -1.071263908506e-6f,1.83128394654e-7f,-3.1353509361e-8f,
	    5.372808776e-9f,-9.21168141e-10f,1.57981265e-10f,-2.7098646e-11f,
	    4.648722e-12f,-7.97527e-13f,1.36827e-13f,-2.3475e-14f,4.027e-15f,
	    -6.91e-16f,1.18e-16f,-2e-17f };

__device__ __constant__ real apsics[16] = { -.0204749044678185f,-.0101801271534859f,
	    5.59718725387e-5f,-1.291717657e-6f,5.72858606e-8f,-3.8213539e-9f,
	    3.397434e-10f,-3.74838e-11f,4.899e-12f,-7.344e-13f,1.233e-13f,
	    -2.28e-14f,4.5e-15f,-9e-16f,2e-16f,-0.f };

/* DECK PSI */
__device__ real slatec_psi(real *x)
{
    /* Initialized data */

    const real pi = 3.14159265358979324f;
     logical first = TRUE_;

    /* System generated locals */
    integer i__1;
    real ret_val, r__1, r__2;

    /* Builtin functions */
    __device__ double sqrt(doublereal), r_int(real *), log(doublereal);

    /* Local variables */
     integer i__, n;
     real y;
    __device__ extern doublereal cot_(real *);
     real aux, xbig;
    __device__ extern doublereal csevl_(real *, real *, const integer *);
     real dxrel;
    __device__ extern integer inits_(real *, const integer *, real *);
     integer ntpsi;
    __device__ extern doublereal r1mach_(const integer *);
     integer ntapsi;
    __device__ extern /* Subroutine */ int xermsg_(char *, char *, char *, const integer *, 
	    const integer *, ftnlen, ftnlen, ftnlen);

/* ***BEGIN PROLOGUE  PSI */
/* ***PURPOSE  Compute the Psi (or Digamma) function. */
/* ***LIBRARY   SLATEC (FNLIB) */
/* ***CATEGORY  C7C */
/* ***TYPE      SINGLE PRECISION (PSI-S, DPSI-D, CPSI-C) */
/* ***KEYWORDS  DIGAMMA FUNCTION, FNLIB, PSI FUNCTION, SPECIAL FUNCTIONS */
/* ***AUTHOR  Fullerton, W., (LANL) */
/* ***DESCRIPTION */

/* PSI(X) calculates the psi (or digamma) function for real argument X. */
/* PSI(X) is the logarithmic derivative of the gamma function of X. */

/* Series for PSI        on the interval  0.          to  1.00000D+00 */
/*                                        with weighted error   2.03E-17 */
/*                                         log weighted error  16.69 */
/*                               significant figures required  16.39 */
/*                                    decimal places required  17.37 */

/* Series for APSI       on the interval  0.          to  2.50000D-01 */
/*                                        with weighted error   5.54E-17 */
/*                                         log weighted error  16.26 */
/*                               significant figures required  14.42 */
/*                                    decimal places required  16.86 */

/* ***REFERENCES  (NONE) */
/* ***ROUTINES CALLED  COT, CSEVL, INITS, R1MACH, XERMSG */
/* ***REVISION HISTORY  (YYMMDD) */
/*   770401  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890531  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ) */
/*   900727  Added EXTERNAL statement.  (WRB) */
/*   920618  Removed space from variable names.  (RWC, WRB) */
/* ***END PROLOGUE  PSI */
/* ***FIRST EXECUTABLE STATEMENT  PSI */
    if (first) {
	r__1 = r1mach_(&c__3) * .1f;
	ntpsi = inits_(psics, &c__23, &r__1);
	r__1 = r1mach_(&c__3) * .1f;
	ntapsi = inits_(apsics, &c__16, &r__1);

	xbig = 1.f / sqrt(r1mach_(&c__3));
	dxrel = sqrt(r1mach_(&c__4));
    }
    first = FALSE_;

    y = dabs(*x);
    if (y >= 2.f) {
	goto L30;
    }

/* PSI(X) FOR -2. .LT. X .LT. 2. */

    n = *x;
    if (*x < 0.f) {
	--n;
    }
    y = *x - n;
    --n;
    r__1 = y * 2.f - 1.f;
    ret_val = csevl_(&r__1, psics, &ntpsi);
    if (n == 0) {
	return ret_val;
    }

    n = -n;
    if (*x == 0.f) {
	xermsg_("SLATEC", "PSI", "X IS 0", &c__2, &c__2, (ftnlen)6, (ftnlen)3,
		 (ftnlen)6);
    }
    if (*x < 0.f && *x + n - 2 == 0.f) {
	xermsg_("SLATEC", "PSI", "X IS A NEGATIVE INTEGER", &c__3, &c__2, (
		ftnlen)6, (ftnlen)3, (ftnlen)23);
    }
    r__2 = *x - .5f;
    if (*x < -.5f && (r__1 = (*x - r_int(&r__2)) / *x, dabs(r__1)) < dxrel) {
	xermsg_("SLATEC", "PSI", "ANSWER LT HALF PRECISION BECAUSE X TOO NEA"
		"R NEGATIVE INTEGER", &c__1, &c__1, (ftnlen)6, (ftnlen)3, (
		ftnlen)60);
    }

    i__1 = n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	ret_val -= 1.f / (*x + i__ - 1);
/* L20: */
    }
    return ret_val;

/* PSI(X) FOR ABS(X) .GE. 2. */

L30:
    aux = 0.f;
    if (y < xbig) {
/* Computing 2nd power */
	r__2 = y;
	r__1 = 8.f / (r__2 * r__2) - 1.f;
	aux = csevl_(&r__1, apsics, &ntapsi);
    }
    if (*x < 0.f) {
	r__1 = pi * *x;
	ret_val = log((dabs(*x))) - .5f / *x + aux - pi * cot_(&r__1);
    }
    if (*x > 0.f) {
	ret_val = log(*x) - .5f / *x + aux;
    }
    return ret_val;

} /* psi_ */

