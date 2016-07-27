/* dpsifn.f -- translated by f2c (version 20100827).
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

static integer c__15 = 15;
static integer c__16 = 16;
static integer c__5 = 5;
static integer c__4 = 4;
static integer c__14 = 14;

/* DECK DPSIFN */
/* Subroutine */ int dpsifn_(doublereal *x, integer *n, integer *kode, 
	integer *m, doublereal *ans, integer *nz, integer *ierr)
{
    /* Initialized data */

    static integer nmax = 100;
    static doublereal b[22] = { 1.,-.5,.166666666666666667,
	    -.0333333333333333333,.0238095238095238095,-.0333333333333333333,
	    .0757575757575757576,-.253113553113553114,1.16666666666666667,
	    -7.09215686274509804,54.9711779448621554,-529.124242424242424,
	    6192.1231884057971,-86580.2531135531136,1425517.16666666667,
	    -27298231.067816092,601580873.900642368,-15116315767.0921569,
	    429614643061.166667,-13711655205088.3328,488332318973593.167,
	    -19296579341940068.1 };

    /* System generated locals */
    integer i__1, i__2;
    doublereal d__1, d__2;

    /* Builtin functions */
    double log(doublereal), exp(doublereal), pow_di(doublereal *, integer *);

    /* Local variables */
    static integer i__, j, k;
    static doublereal s, t, t1, t2;
    static integer fn;
    static doublereal ta;
    static integer mm, nn, np;
    static doublereal fx, tk;
    static integer mx, nx;
    static doublereal xm, tt, xq, den, arg, fln, r1m4, r1m5, eps, rln, tol, 
	    xln, trm[22], tss, tst, elim, xinc, xmin, tols, xdmy, yint, trmr[
	    100], rxsq, slope, xdmln, wdtol;
    extern doublereal d1mach_(integer *);
    extern integer i1mach_(integer *);

/* ***BEGIN PROLOGUE  DPSIFN */
/* ***PURPOSE  Compute derivatives of the Psi function. */
/* ***LIBRARY   SLATEC */
/* ***CATEGORY  C7C */
/* ***TYPE      DOUBLE PRECISION (PSIFN-S, DPSIFN-D) */
/* ***KEYWORDS  DERIVATIVES OF THE GAMMA FUNCTION, POLYGAMMA FUNCTION, */
/*             PSI FUNCTION */
/* ***AUTHOR  Amos, D. E., (SNLA) */
/* ***DESCRIPTION */

/*         The following definitions are used in DPSIFN: */

/*      Definition 1 */
/*         PSI(X) = d/dx (ln(GAMMA(X)), the first derivative of */
/*                  the log GAMMA function. */
/*      Definition 2 */
/*                     K   K */
/*         PSI(K,X) = d /dx (PSI(X)), the K-th derivative of PSI(X). */
/*   ___________________________________________________________________ */
/*      DPSIFN computes a sequence of SCALED derivatives of */
/*      the PSI function; i.e. for fixed X and M it computes */
/*      the M-member sequence */

/*                    ((-1)**(K+1)/GAMMA(K+1))*PSI(K,X) */
/*                       for K = N,...,N+M-1 */

/*      where PSI(K,X) is as defined above.   For KODE=1, DPSIFN returns */
/*      the scaled derivatives as described.  KODE=2 is operative only */
/*      when K=0 and in that case DPSIFN returns -PSI(X) + LN(X).  That */
/*      is, the logarithmic behavior for large X is removed when KODE=2 */
/*      and K=0.  When sums or differences of PSI functions are computed */
/*      the logarithmic terms can be combined analytically and computed */
/*      separately to help retain significant digits. */

/*         Note that CALL DPSIFN(X,0,1,1,ANS) results in */
/*                   ANS = -PSI(X) */

/*     Input      X is DOUBLE PRECISION */
/*           X      - Argument, X .gt. 0.0D0 */
/*           N      - First member of the sequence, 0 .le. N .le. 100 */
/*                    N=0 gives ANS(1) = -PSI(X)       for KODE=1 */
/*                                       -PSI(X)+LN(X) for KODE=2 */
/*           KODE   - Selection parameter */
/*                    KODE=1 returns scaled derivatives of the PSI */
/*                    function. */
/*                    KODE=2 returns scaled derivatives of the PSI */
/*                    function EXCEPT when N=0. In this case, */
/*                    ANS(1) = -PSI(X) + LN(X) is returned. */
/*           M      - Number of members of the sequence, M.ge.1 */

/*    Output     ANS is DOUBLE PRECISION */
/*           ANS    - A vector of length at least M whose first M */
/*                    components contain the sequence of derivatives */
/*                    scaled according to KODE. */
/*           NZ     - Underflow flag */
/*                    NZ.eq.0, A normal return */
/*                    NZ.ne.0, Underflow, last NZ components of ANS are */
/*                             set to zero, ANS(M-K+1)=0.0, K=1,...,NZ */
/*           IERR   - Error flag */
/*                    IERR=0, A normal return, computation completed */
/*                    IERR=1, Input error,     no computation */
/*                    IERR=2, Overflow,        X too small or N+M-1 too */
/*                            large or both */
/*                    IERR=3, Error,           N too large. Dimensioned */
/*                            array TRMR(NMAX) is not large enough for N */

/*         The nominal computational accuracy is the maximum of unit */
/*         roundoff (=D1MACH(4)) and 1.0D-18 since critical constants */
/*         are given to only 18 digits. */

/*         PSIFN is the single precision version of DPSIFN. */

/* *Long Description: */

/*         The basic method of evaluation is the asymptotic expansion */
/*         for large X.ge.XMIN followed by backward recursion on a two */
/*         term recursion relation */

/*                  W(X+1) + X**(-N-1) = W(X). */

/*         This is supplemented by a series */

/*                  SUM( (X+K)**(-N-1) , K=0,1,2,... ) */

/*         which converges rapidly for large N. Both XMIN and the */
/*         number of terms of the series are calculated from the unit */
/*         roundoff of the machine environment. */

/* ***REFERENCES  Handbook of Mathematical Functions, National Bureau */
/*                 of Standards Applied Mathematics Series 55, edited */
/*                 by M. Abramowitz and I. A. Stegun, equations 6.3.5, */
/*                 6.3.18, 6.4.6, 6.4.9 and 6.4.10, pp.258-260, 1964. */
/*               D. E. Amos, A portable Fortran subroutine for */
/*                 derivatives of the Psi function, Algorithm 610, ACM */
/*                 Transactions on Mathematical Software 9, 4 (1983), */
/*                 pp. 494-502. */
/* ***ROUTINES CALLED  D1MACH, I1MACH */
/* ***REVISION HISTORY  (YYMMDD) */
/*   820601  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890911  Removed unnecessary intrinsics.  (WRB) */
/*   891006  Cosmetic changes to prologue.  (WRB) */
/*   891006  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   920501  Reformatted the REFERENCES section.  (WRB) */
/* ***END PROLOGUE  DPSIFN */
    /* Parameter adjustments */
    --ans;

    /* Function Body */
/* ----------------------------------------------------------------------- */
/*             BERNOULLI NUMBERS */
/* ----------------------------------------------------------------------- */

/* ***FIRST EXECUTABLE STATEMENT  DPSIFN */
    *ierr = 0;
    *nz = 0;
    if (*x <= 0.) {
	*ierr = 1;
    }
    if (*n < 0) {
	*ierr = 1;
    }
    if (*kode < 1 || *kode > 2) {
	*ierr = 1;
    }
    if (*m < 1) {
	*ierr = 1;
    }
    if (*ierr != 0) {
	return 0;
    }
    mm = *m;
/* Computing MIN */
    i__1 = -i1mach_(&c__15), i__2 = i1mach_(&c__16);
    nx = min(i__1,i__2);
    r1m5 = d1mach_(&c__5);
    r1m4 = d1mach_(&c__4) * .5;
    wdtol = max(r1m4,5e-19);
/* ----------------------------------------------------------------------- */
/*     ELIM = APPROXIMATE EXPONENTIAL OVER AND UNDERFLOW LIMIT */
/* ----------------------------------------------------------------------- */
    elim = (nx * r1m5 - 3.) * 2.302;
    xln = log(*x);
L41:
    nn = *n + mm - 1;
    fn = nn;
    t = (fn + 1) * xln;
/* ----------------------------------------------------------------------- */
/*     OVERFLOW AND UNDERFLOW TEST FOR SMALL AND LARGE X */
/* ----------------------------------------------------------------------- */
    if (abs(t) > elim) {
	goto L290;
    }
    if (*x < wdtol) {
	goto L260;
    }
/* ----------------------------------------------------------------------- */
/*     COMPUTE XMIN AND THE NUMBER OF TERMS OF THE SERIES, FLN+1 */
/* ----------------------------------------------------------------------- */
    rln = r1m5 * i1mach_(&c__14);
    rln = min(rln,18.06);
    fln = max(rln,3.) - 3.;
    yint = fln * .4 + 3.5;
    slope = fln * (fln * 6.038e-4 + .008677) + .21;
    xm = yint + slope * fn;
    mx = (integer) xm + 1;
    xmin = (doublereal) mx;
    if (*n == 0) {
	goto L50;
    }
    xm = rln * -2.302 - min(0.,xln);
    arg = xm / *n;
    arg = min(0.,arg);
    eps = exp(arg);
    xm = 1. - eps;
    if (abs(arg) < .001) {
	xm = -arg;
    }
    fln = *x * xm / eps;
    xm = xmin - *x;
    if (xm > 7. && fln < 15.) {
	goto L200;
    }
L50:
    xdmy = *x;
    xdmln = xln;
    xinc = 0.;
    if (*x >= xmin) {
	goto L60;
    }
    nx = (integer) (*x);
    xinc = xmin - nx;
    xdmy = *x + xinc;
    xdmln = log(xdmy);
L60:
/* ----------------------------------------------------------------------- */
/*     GENERATE W(N+MM-1,X) BY THE ASYMPTOTIC EXPANSION */
/* ----------------------------------------------------------------------- */
    t = fn * xdmln;
    t1 = xdmln + xdmln;
    t2 = t + xdmln;
/* Computing MAX */
    d__1 = abs(t), d__2 = abs(t1), d__1 = max(d__1,d__2), d__2 = abs(t2);
    tk = max(d__1,d__2);
    if (tk > elim) {
	goto L380;
    }
    tss = exp(-t);
    tt = .5 / xdmy;
    t1 = tt;
    tst = wdtol * tt;
    if (nn != 0) {
	t1 = tt + 1. / fn;
    }
    rxsq = 1. / (xdmy * xdmy);
    ta = rxsq * .5;
    t = (fn + 1) * ta;
    s = t * b[2];
    if (abs(s) < tst) {
	goto L80;
    }
    tk = 2.;
    for (k = 4; k <= 22; ++k) {
	t = t * ((tk + fn + 1) / (tk + 1.)) * ((tk + fn) / (tk + 2.)) * rxsq;
	trm[k - 1] = t * b[k - 1];
	if ((d__1 = trm[k - 1], abs(d__1)) < tst) {
	    goto L80;
	}
	s += trm[k - 1];
	tk += 2.;
/* L70: */
    }
L80:
    s = (s + t1) * tss;
    if (xinc == 0.) {
	goto L100;
    }
/* ----------------------------------------------------------------------- */
/*     BACKWARD RECUR FROM XDMY TO X */
/* ----------------------------------------------------------------------- */
    nx = (integer) xinc;
    np = nn + 1;
    if (nx > nmax) {
	goto L390;
    }
    if (nn == 0) {
	goto L160;
    }
    xm = xinc - 1.;
    fx = *x + xm;
/* ----------------------------------------------------------------------- */
/*     THIS LOOP SHOULD NOT BE CHANGED. FX IS ACCURATE WHEN X IS SMALL */
/* ----------------------------------------------------------------------- */
    i__1 = nx;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = -np;
	trmr[i__ - 1] = pow_di(&fx, &i__2);
	s += trmr[i__ - 1];
	xm += -1.;
	fx = *x + xm;
/* L90: */
    }
L100:
    ans[mm] = s;
    if (fn == 0) {
	goto L180;
    }
/* ----------------------------------------------------------------------- */
/*     GENERATE LOWER DERIVATIVES, J.LT.N+MM-1 */
/* ----------------------------------------------------------------------- */
    if (mm == 1) {
	return 0;
    }
    i__1 = mm;
    for (j = 2; j <= i__1; ++j) {
	--fn;
	tss *= xdmy;
	t1 = tt;
	if (fn != 0) {
	    t1 = tt + 1. / fn;
	}
	t = (fn + 1) * ta;
	s = t * b[2];
	if (abs(s) < tst) {
	    goto L120;
	}
	tk = (doublereal) (fn + 4);
	for (k = 4; k <= 22; ++k) {
	    trm[k - 1] = trm[k - 1] * (fn + 1) / tk;
	    if ((d__1 = trm[k - 1], abs(d__1)) < tst) {
		goto L120;
	    }
	    s += trm[k - 1];
	    tk += 2.;
/* L110: */
	}
L120:
	s = (s + t1) * tss;
	if (xinc == 0.) {
	    goto L140;
	}
	if (fn == 0) {
	    goto L160;
	}
	xm = xinc - 1.;
	fx = *x + xm;
	i__2 = nx;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    trmr[i__ - 1] *= fx;
	    s += trmr[i__ - 1];
	    xm += -1.;
	    fx = *x + xm;
/* L130: */
	}
L140:
	mx = mm - j + 1;
	ans[mx] = s;
	if (fn == 0) {
	    goto L180;
	}
/* L150: */
    }
    return 0;
/* ----------------------------------------------------------------------- */
/*     RECURSION FOR N = 0 */
/* ----------------------------------------------------------------------- */
L160:
    i__1 = nx;
    for (i__ = 1; i__ <= i__1; ++i__) {
	s += 1. / (*x + nx - i__);
/* L170: */
    }
L180:
    if (*kode == 2) {
	goto L190;
    }
    ans[1] = s - xdmln;
    return 0;
L190:
    if (xdmy == *x) {
	return 0;
    }
    xq = xdmy / *x;
    ans[1] = s - log(xq);
    return 0;
/* ----------------------------------------------------------------------- */
/*     COMPUTE BY SERIES (X+K)**(-(N+1)) , K=0,1,2,... */
/* ----------------------------------------------------------------------- */
L200:
    nn = (integer) fln + 1;
    np = *n + 1;
    t1 = (*n + 1) * xln;
    t = exp(-t1);
    s = t;
    den = *x;
    i__1 = nn;
    for (i__ = 1; i__ <= i__1; ++i__) {
	den += 1.;
	i__2 = -np;
	trm[i__ - 1] = pow_di(&den, &i__2);
	s += trm[i__ - 1];
/* L210: */
    }
    ans[1] = s;
    if (*n != 0) {
	goto L220;
    }
    if (*kode == 2) {
	ans[1] = s + xln;
    }
L220:
    if (mm == 1) {
	return 0;
    }
/* ----------------------------------------------------------------------- */
/*     GENERATE HIGHER DERIVATIVES, J.GT.N */
/* ----------------------------------------------------------------------- */
    tol = wdtol / 5.;
    i__1 = mm;
    for (j = 2; j <= i__1; ++j) {
	t /= *x;
	s = t;
	tols = t * tol;
	den = *x;
	i__2 = nn;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    den += 1.;
	    trm[i__ - 1] /= den;
	    s += trm[i__ - 1];
	    if (trm[i__ - 1] < tols) {
		goto L240;
	    }
/* L230: */
	}
L240:
	ans[j] = s;
/* L250: */
    }
    return 0;
/* ----------------------------------------------------------------------- */
/*     SMALL X.LT.UNIT ROUND OFF */
/* ----------------------------------------------------------------------- */
L260:
    i__1 = -(*n) - 1;
    ans[1] = pow_di(x, &i__1);
    if (mm == 1) {
	goto L280;
    }
    k = 1;
    i__1 = mm;
    for (i__ = 2; i__ <= i__1; ++i__) {
	ans[k + 1] = ans[k] / *x;
	++k;
/* L270: */
    }
L280:
    if (*n != 0) {
	return 0;
    }
    if (*kode == 2) {
	ans[1] += xln;
    }
    return 0;
L290:
    if (t > 0.) {
	goto L380;
    }
    *nz = 0;
    *ierr = 2;
    return 0;
L380:
    ++(*nz);
    ans[mm] = 0.;
    --mm;
    if (mm == 0) {
	return 0;
    }
    goto L41;
L390:
    *nz = 0;
    *ierr = 3;
    return 0;
} /* dpsifn_ */

