/* cot.f -- translated by f2c (version 20100827).
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

static integer c__3 = 3;
static integer c__8 = 8;
static integer c__4 = 4;
static integer c__2 = 2;
static integer c__1 = 1;
static real c_b19 = 2.f;

/* DECK COT */
doublereal cot_(real *x)
{
    /* Initialized data */

    static real cotcs[8] = { .2402591609829563f,-.016533031601500228f,
	    -4.2998391931724e-5f,-1.59283223327e-7f,-6.19109313e-10f,
	    -2.430197e-12f,-9.56e-15f,-3.7e-17f };
    static real pi2rec = .011619772367581343f;
    static logical first = TRUE_;

    /* System generated locals */
    real ret_val, r__1, r__2;

    /* Builtin functions */
    double sqrt(doublereal), log(doublereal), exp(doublereal), r_int(real *), 
	    r_mod(real *, real *), r_sign(real *, real *);

    /* Local variables */
    static real y;
    static integer ifn;
    static real xmin, yrem, xmax, xsml;
    extern doublereal csevl_(real *, real *, integer *);
    static real ainty;
    extern integer inits_(real *, integer *, real *);
    static real sqeps;
    extern doublereal r1mach_(integer *);
    static real ainty2, prodbg;
    extern /* Subroutine */ int xermsg_(char *, char *, char *, integer *, 
	    integer *, ftnlen, ftnlen, ftnlen);
    static integer nterms;

/* ***BEGIN PROLOGUE  COT */
/* ***PURPOSE  Compute the cotangent. */
/* ***LIBRARY   SLATEC (FNLIB) */
/* ***CATEGORY  C4A */
/* ***TYPE      SINGLE PRECISION (COT-S, DCOT-D, CCOT-C) */
/* ***KEYWORDS  COTANGENT, ELEMENTARY FUNCTIONS, FNLIB, TRIGONOMETRIC */
/* ***AUTHOR  Fullerton, W., (LANL) */
/* ***DESCRIPTION */

/* COT(X) calculates the cotangent of the real argument X.  X is in */
/* units of radians. */

/* Series for COT        on the interval  0.          to  6.25000D-02 */
/*                                        with weighted error   3.76E-17 */
/*                                         log weighted error  16.42 */
/*                               significant figures required  15.51 */
/*                                    decimal places required  16.88 */

/* ***REFERENCES  (NONE) */
/* ***ROUTINES CALLED  CSEVL, INITS, R1MACH, XERMSG */
/* ***REVISION HISTORY  (YYMMDD) */
/*   770601  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890531  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ) */
/*   920618  Removed space from variable names.  (RWC, WRB) */
/* ***END PROLOGUE  COT */
/* ***FIRST EXECUTABLE STATEMENT  COT */
    if (first) {
	r__1 = r1mach_(&c__3) * .1f;
	nterms = inits_(cotcs, &c__8, &r__1);
	xmax = 1.f / r1mach_(&c__4);
	xsml = sqrt(r1mach_(&c__3) * 3.f);
/* Computing MAX */
	r__1 = log(r1mach_(&c__1)), r__2 = -log(r1mach_(&c__2));
	xmin = exp(dmax(r__1,r__2) + .01f);
	sqeps = sqrt(r1mach_(&c__4));
    }
    first = FALSE_;

    y = dabs(*x);
    if (dabs(*x) < xmin) {
	xermsg_("SLATEC", "COT", "ABS(X) IS ZERO OR SO SMALL COT OVERFLOWS", &
		c__2, &c__2, (ftnlen)6, (ftnlen)3, (ftnlen)40);
    }
    if (y > xmax) {
	xermsg_("SLATEC", "COT", "NO PRECISION BECAUSE ABS(X) IS TOO BIG", &
		c__3, &c__2, (ftnlen)6, (ftnlen)3, (ftnlen)38);
    }

/* CAREFULLY COMPUTE Y * (2/PI) = (AINT(Y) + REM(Y)) * (.625 + PI2REC) */
/* = AINT(.625*Y) + REM(.625*Y) + Y*PI2REC  =  AINT(.625*Y) + Z */
/* = AINT(.625*Y) + AINT(Z) + REM(Z) */

    ainty = r_int(&y);
    yrem = y - ainty;
    prodbg = ainty * .625f;
    ainty = r_int(&prodbg);
    y = prodbg - ainty + yrem * .625f + y * pi2rec;
    ainty2 = r_int(&y);
    ainty += ainty2;
    y -= ainty2;

    ifn = r_mod(&ainty, &c_b19);
    if (ifn == 1) {
	y = 1.f - y;
    }

    if (dabs(*x) > .5f && y < dabs(*x) * sqeps) {
	xermsg_("SLATEC", "COT", "ANSWER LT HALF PRECISION, ABS(X) TOO BIG O"
		"R X NEAR N*PI (N.NE.0)", &c__1, &c__1, (ftnlen)6, (ftnlen)3, (
		ftnlen)64);
    }

    if (y > .25f) {
	goto L20;
    }
    ret_val = 1.f / *x;
    if (y > xsml) {
	r__1 = y * 32.f * y - 1.f;
	ret_val = (csevl_(&r__1, cotcs, &nterms) + .5f) / y;
    }
    goto L40;

L20:
    if (y > .5f) {
	goto L30;
    }
    r__1 = y * 8.f * y - 1.f;
    ret_val = (csevl_(&r__1, cotcs, &nterms) + .5f) / (y * .5f);
/* Computing 2nd power */
    r__1 = ret_val;
    ret_val = (r__1 * r__1 - 1.f) * .5f / ret_val;
    goto L40;

L30:
    r__1 = y * 2.f * y - 1.f;
    ret_val = (csevl_(&r__1, cotcs, &nterms) + .5f) / (y * .25f);
/* Computing 2nd power */
    r__1 = ret_val;
    ret_val = (r__1 * r__1 - 1.f) * .5f / ret_val;
/* Computing 2nd power */
    r__1 = ret_val;
    ret_val = (r__1 * r__1 - 1.f) * .5f / ret_val;

L40:
    if (*x != 0.f) {
	ret_val = r_sign(&ret_val, x);
    }
    if (ifn == 1) {
	ret_val = -ret_val;
    }

    return ret_val;
} /* cot_ */

