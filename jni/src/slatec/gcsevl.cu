/* csevl.f -- translated by f2c (version 20100827).
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

/* DECK CSEVL */
__device__ __forceinline__ doublereal csevl_(real *x, real *cs, const integer *n)
{
    /* Initialized data */

     logical first = TRUE_;

    /* System generated locals */
    integer i__1;
    real ret_val;

    /* Local variables */
     integer i__;
     real b0, b1, b2;
     integer ni;
     real twox, onepl;

/* ***BEGIN PROLOGUE  CSEVL */
/* ***PURPOSE  Evaluate a Chebyshev series. */
/* ***LIBRARY   SLATEC (FNLIB) */
/* ***CATEGORY  C3A2 */
/* ***TYPE      SINGLE PRECISION (CSEVL-S, DCSEVL-D) */
/* ***KEYWORDS  CHEBYSHEV SERIES, FNLIB, SPECIAL FUNCTIONS */
/* ***AUTHOR  Fullerton, W., (LANL) */
/* ***DESCRIPTION */

/*  Evaluate the N-term Chebyshev series CS at X.  Adapted from */
/*  a method presented in the paper by Broucke referenced below. */

/*       Input Arguments -- */
/*  X    value at which the series is to be evaluated. */
/*  CS   array of N terms of a Chebyshev series.  In evaluating */
/*       CS, only half the first coefficient is summed. */
/*  N    number of terms in array CS. */

/* ***REFERENCES  R. Broucke, Ten subroutines for the manipulation of */
/*                 Chebyshev series, Algorithm 446, Communications of */
/*                 the A.C.M. 16, (1973) pp. 254-256. */
/*               L. Fox and I. B. Parker, Chebyshev Polynomials in */
/*                 Numerical Analysis, Oxford University Press, 1968, */
/*                 page 56. */
/* ***ROUTINES CALLED  R1MACH, XERMSG */
/* ***REVISION HISTORY  (YYMMDD) */
/*   770401  DATE WRITTEN */
/*   890831  Modified array declarations.  (WRB) */
/*   890831  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ) */
/*   900329  Prologued revised extensively and code rewritten to allow */
/*           X to be slightly outside interval (-1,+1).  (WRB) */
/*   920501  Reformatted the REFERENCES section.  (WRB) */
/* ***END PROLOGUE  CSEVL */
    /* Parameter adjustments */
    --cs;

    /* Function Body */
/* ***FIRST EXECUTABLE STATEMENT  CSEVL */
    if (first) {
	onepl = 1.f + r1mach_(&c__4);
    }
    first = FALSE_;
    if (*n < 1) {
	xermsg_("SLATEC", "CSEVL", "NUMBER OF TERMS .LE. 0", &c__2, &c__2, (
		ftnlen)6, (ftnlen)5, (ftnlen)22);
    }
    if (*n > 1000) {
	xermsg_("SLATEC", "CSEVL", "NUMBER OF TERMS .GT. 1000", &c__3, &c__2, 
		(ftnlen)6, (ftnlen)5, (ftnlen)25);
    }
    if (dabs(*x) > onepl) {
	xermsg_("SLATEC", "CSEVL", "X OUTSIDE THE INTERVAL (-1,+1)", &c__1, &
		c__1, (ftnlen)6, (ftnlen)5, (ftnlen)30);
    }

    b1 = 0.f;
    b0 = 0.f;
    twox = *x * 2.f;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	b2 = b1;
	b1 = b0;
	ni = *n + 1 - i__;
	b0 = twox * b1 - b2 + cs[ni];
/* L10: */
    }

    ret_val = (b0 - b2) * .5f;

    return ret_val;
} /* csevl_ */

