/* dpsi.f -- translated by f2c (version 20100827).
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
static integer c__42 = 42;
static integer c__16 = 16;
static integer c__4 = 4;
static integer c__2 = 2;
static integer c__1 = 1;

/* DECK DPSI */
doublereal dpsi_(doublereal *x)
{
    /* Initialized data */

    static doublereal psics[42] = { -.038057080835217921520437677667039,
	    .49141539302938712748204699654277,
	    -.056815747821244730242892064734081,
	    .0083578212259143131362775650747862,
	    -.0013332328579943425998079274172393,
	    2.2031328706930824892872397979521e-4,
	    -3.7040238178456883592889086949229e-5,
	    6.283793654854989893365141871769e-6,
	    -1.0712639085061849855283541747074e-6,
	    1.8312839465484165805731589810378e-7,
	    -3.1353509361808509869005779796885e-8,
	    5.3728087762007766260471919143615e-9,
	    -9.211681415978427571788063262473e-10,
	    1.5798126521481822782252884032823e-10,
	    -2.7098646132380443065440589409707e-11,
	    4.6487228599096834872947319529549e-12,
	    -7.9752725638303689726504797772737e-13,
	    1.3682723857476992249251053892838e-13,
	    -2.3475156060658972717320677980719e-14,
	    4.0276307155603541107907925006281e-15,
	    -6.9102518531179037846547422974771e-16,
	    1.1856047138863349552929139525768e-16,
	    -2.0341689616261559308154210484223e-17,
	    3.4900749686463043850374232932351e-18,
	    -5.9880146934976711003011081393493e-19,
	    1.0273801628080588258398005712213e-19,
	    -1.7627049424561071368359260105386e-20,
	    3.0243228018156920457454035490133e-21,
	    -5.1889168302092313774286088874666e-22,
	    8.9027730345845713905005887487999e-23,
	    -1.5274742899426728392894971904e-23,
	    2.6207314798962083136358318079999e-24,
	    -4.4964642738220696772598388053333e-25,
	    7.7147129596345107028919364266666e-26,
	    -1.3236354761887702968102638933333e-26,
	    2.2709994362408300091277311999999e-27,
	    -3.8964190215374115954491391999999e-28,
	    6.6851981388855302310679893333333e-29,
	    -1.1469986654920864872529919999999e-29,
	    1.9679385886541405920515413333333e-30,
	    -3.37644881897509798019072e-31,
	    5.7930703193214159246677333333333e-32 };
    static doublereal apsics[16] = { -8.32710791069290760174456932269e-4,
	    -4.1625184219273935282162712199e-4,
	    1.03431560978741291174463193961e-7,
	    -1.21468184135904152987299556365e-10,
	    3.11369431998356155521240278178e-13,
	    -1.36461337193177041776516100945e-15,
	    9.02051751315416565130837974e-18,
	    -8.31542997421591464829933635466e-20,
	    1.01224257073907254188479482666e-21,
	    -1.56270249435622507620478933333e-23,
	    2.96542716808903896133226666666e-25,
	    -6.74686886765702163741866666666e-27,
	    1.80345311697189904213333333333e-28,
	    -5.56901618245983607466666666666e-30,
	    1.95867922607736251733333333333e-31,-7.7519589252333568e-33 };
    static doublereal pi = 3.1415926535897932384626433832795;
    static logical first = TRUE_;

    /* System generated locals */
    integer i__1;
    real r__1;
    doublereal ret_val, d__1, d__2;

    /* Builtin functions */
    double sqrt(doublereal), d_int(doublereal *), log(doublereal);

    /* Local variables */
    static integer i__, n;
    static doublereal y, aux, xbig;
    extern doublereal dcot_(doublereal *);
    static doublereal dxrel;
    extern doublereal d1mach_(integer *);
    static integer ntpsi;
    extern doublereal dcsevl_(doublereal *, doublereal *, integer *);
    extern integer initds_(doublereal *, integer *, real *);
    static integer ntapsi;
    extern /* Subroutine */ int xermsg_(char *, char *, char *, integer *, 
	    integer *, ftnlen, ftnlen, ftnlen);

/* ***BEGIN PROLOGUE  DPSI */
/* ***PURPOSE  Compute the Psi (or Digamma) function. */
/* ***LIBRARY   SLATEC (FNLIB) */
/* ***CATEGORY  C7C */
/* ***TYPE      DOUBLE PRECISION (PSI-S, DPSI-D, CPSI-C) */
/* ***KEYWORDS  DIGAMMA FUNCTION, FNLIB, PSI FUNCTION, SPECIAL FUNCTIONS */
/* ***AUTHOR  Fullerton, W., (LANL) */
/* ***DESCRIPTION */

/* DPSI calculates the double precision Psi (or Digamma) function for */
/* double precision argument X.  PSI(X) is the logarithmic derivative */
/* of the Gamma function of X. */

/* Series for PSI        on the interval  0.          to  1.00000E+00 */
/*                                        with weighted error   5.79E-32 */
/*                                         log weighted error  31.24 */
/*                               significant figures required  30.93 */
/*                                    decimal places required  32.05 */


/* Series for APSI       on the interval  0.          to  1.00000E-02 */
/*                                        with weighted error   7.75E-33 */
/*                                         log weighted error  32.11 */
/*                               significant figures required  28.88 */
/*                                    decimal places required  32.71 */

/* ***REFERENCES  (NONE) */
/* ***ROUTINES CALLED  D1MACH, DCOT, DCSEVL, INITDS, XERMSG */
/* ***REVISION HISTORY  (YYMMDD) */
/*   770601  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890911  Removed unnecessary intrinsics.  (WRB) */
/*   890911  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ) */
/*   900727  Added EXTERNAL statement.  (WRB) */
/*   920618  Removed space from variable name.  (RWC, WRB) */
/* ***END PROLOGUE  DPSI */
/* ***FIRST EXECUTABLE STATEMENT  DPSI */
    if (first) {
	r__1 = (real) d1mach_(&c__3) * .1f;
	ntpsi = initds_(psics, &c__42, &r__1);
	r__1 = (real) d1mach_(&c__3) * .1f;
	ntapsi = initds_(apsics, &c__16, &r__1);

	xbig = 1. / sqrt(d1mach_(&c__3));
	dxrel = sqrt(d1mach_(&c__4));
    }
    first = FALSE_;

    y = abs(*x);

    if (y > 10.) {
	goto L50;
    }

/* DPSI(X) FOR ABS(X) .LE. 2 */

    n = (integer) (*x);
    if (*x < 0.) {
	--n;
    }
    y = *x - n;
    --n;
    d__1 = y * 2. - 1.;
    ret_val = dcsevl_(&d__1, psics, &ntpsi);
    if (n == 0) {
	return ret_val;
    }

    if (n > 0) {
	goto L30;
    }

    n = -n;
    if (*x == 0.) {
	xermsg_("SLATEC", "DPSI", "X IS 0", &c__2, &c__2, (ftnlen)6, (ftnlen)
		4, (ftnlen)6);
    }
    if (*x < 0. && *x + n - 2 == 0.) {
	xermsg_("SLATEC", "DPSI", "X IS A NEGATIVE INTEGER", &c__3, &c__2, (
		ftnlen)6, (ftnlen)4, (ftnlen)23);
    }
    d__2 = *x - .5;
    if (*x < -.5 && (d__1 = (*x - d_int(&d__2)) / *x, abs(d__1)) < dxrel) {
	xermsg_("SLATEC", "DPSI", "ANSWER LT HALF PRECISION BECAUSE X TOO NE"
		"AR NEGATIVE INTEGER", &c__1, &c__1, (ftnlen)6, (ftnlen)4, (
		ftnlen)60);
    }

    i__1 = n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	ret_val -= 1. / (*x + i__ - 1);
/* L20: */
    }
    return ret_val;

/* DPSI(X) FOR X .GE. 2.0 AND X .LE. 10.0 */

L30:
    i__1 = n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	ret_val += 1. / (y + i__);
/* L40: */
    }
    return ret_val;

/* DPSI(X) FOR ABS(X) .GT. 10.0 */

L50:
    aux = 0.;
    if (y < xbig) {
/* Computing 2nd power */
	d__2 = 10. / y;
	d__1 = d__2 * d__2 * 2. - 1.;
	aux = dcsevl_(&d__1, apsics, &ntapsi);
    }

    if (*x < 0.) {
	d__1 = pi * *x;
	ret_val = log((abs(*x))) - .5 / *x + aux - pi * dcot_(&d__1);
    }
    if (*x > 0.) {
	ret_val = log(*x) - .5 / *x + aux;
    }
    return ret_val;

} /* dpsi_ */

