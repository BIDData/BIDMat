/* r1mach.f -- translated by f2c (version 20100827).
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

static integer c__1 = 1;
static integer c__2 = 2;

const real equiv_4[5] = {
  1.18E-38,
  3.40E+38,
  0.595E-07,
  1.19E-07,
  0.30102999566
};


/* DECK R1MACH */
doublereal r1mach_(integer *i__)
{
    /* System generated locals */
    real ret_val;

    /* Local variables */
#define log10 ((integer *)equiv_4 + 4)
#define large ((integer *)equiv_4 + 1)
#define rmach (equiv_4)
#define small ((integer *)equiv_4)
#define diver ((integer *)equiv_4 + 3)
#define right ((integer *)equiv_4 + 2)
    extern /* Subroutine */ int xermsg_(char *, char *, char *, integer *, 
	    integer *, ftnlen, ftnlen, ftnlen);

/* ***BEGIN PROLOGUE  R1MACH */
/* ***PURPOSE  Return floating point machine dependent constants. */
/* ***LIBRARY   SLATEC */
/* ***CATEGORY  R1 */
/* ***TYPE      SINGLE PRECISION (R1MACH-S, D1MACH-D) */
/* ***KEYWORDS  MACHINE CONSTANTS */
/* ***AUTHOR  Fox, P. A., (Bell Labs) */
/*           Hall, A. D., (Bell Labs) */
/*           Schryer, N. L., (Bell Labs) */
/* ***DESCRIPTION */

/*   R1MACH can be used to obtain machine-dependent parameters for the */
/*   local machine environment.  It is a function subprogram with one */
/*   (input) argument, and can be referenced as follows: */

/*        A = R1MACH(I) */

/*   where I=1,...,5.  The (output) value of A above is determined by */
/*   the (input) value of I.  The results for various values of I are */
/*   discussed below. */

/*   R1MACH(1) = B**(EMIN-1), the smallest positive magnitude. */
/*   R1MACH(2) = B**EMAX*(1 - B**(-T)), the largest magnitude. */
/*   R1MACH(3) = B**(-T), the smallest relative spacing. */
/*   R1MACH(4) = B**(1-T), the largest relative spacing. */
/*   R1MACH(5) = LOG10(B) */

/*   Assume single precision numbers are represented in the T-digit, */
/*   base-B form */

/*              sign (B**E)*( (X(1)/B) + ... + (X(T)/B**T) ) */

/*   where 0 .LE. X(I) .LT. B for I=1,...,T, 0 .LT. X(1), and */
/*   EMIN .LE. E .LE. EMAX. */

/*   The values of B, T, EMIN and EMAX are provided in I1MACH as */
/*   follows: */
/*   I1MACH(10) = B, the base. */
/*   I1MACH(11) = T, the number of base-B digits. */
/*   I1MACH(12) = EMIN, the smallest exponent E. */
/*   I1MACH(13) = EMAX, the largest exponent E. */

/*   To alter this function for a particular environment, the desired */
/*   set of DATA statements should be activated by removing the C from */
/*   column 1.  Also, the values of R1MACH(1) - R1MACH(4) should be */
/*   checked for consistency with the local operating system. */

/* ***REFERENCES  P. A. Fox, A. D. Hall and N. L. Schryer, Framework for */
/*                 a portable library, ACM Transactions on Mathematical */
/*                 Software 4, 2 (June 1978), pp. 177-188. */
/* ***ROUTINES CALLED  XERMSG */
/* ***REVISION HISTORY  (YYMMDD) */
/*   790101  DATE WRITTEN */
/*   890213  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ) */
/*   900618  Added DEC RISC constants.  (WRB) */
/*   900723  Added IBM RS 6000 constants.  (WRB) */
/*   910710  Added HP 730 constants.  (SMR) */
/*   911114  Added Convex IEEE constants.  (WRB) */
/*   920121  Added SUN -r8 compiler option constants.  (WRB) */
/*   920229  Added Touchstone Delta i860 constants.  (WRB) */
/*   920501  Reformatted the REFERENCES section.  (WRB) */
/*   920625  Added CONVEX -p8 and -pd8 compiler option constants. */
/*           (BKS, WRB) */
/*   930201  Added DEC Alpha and SGI constants.  (RWC and WRB) */
/* ***END PROLOGUE  R1MACH */




/*     MACHINE CONSTANTS FOR THE AMIGA */
/*     ABSOFT FORTRAN COMPILER USING THE 68020/68881 COMPILER OPTION */

/*     DATA SMALL(1) / Z'00800000' / */
/*     DATA LARGE(1) / Z'7F7FFFFF' / */
/*     DATA RIGHT(1) / Z'33800000' / */
/*     DATA DIVER(1) / Z'34000000' / */
/*     DATA LOG10(1) / Z'3E9A209B' / */

/*     MACHINE CONSTANTS FOR THE AMIGA */
/*     ABSOFT FORTRAN COMPILER USING SOFTWARE FLOATING POINT */

/*     DATA SMALL(1) / Z'00800000' / */
/*     DATA LARGE(1) / Z'7EFFFFFF' / */
/*     DATA RIGHT(1) / Z'33800000' / */
/*     DATA DIVER(1) / Z'34000000' / */
/*     DATA LOG10(1) / Z'3E9A209B' / */

/*     MACHINE CONSTANTS FOR THE APOLLO */

/*     DATA SMALL(1) / 16#00800000 / */
/*     DATA LARGE(1) / 16#7FFFFFFF / */
/*     DATA RIGHT(1) / 16#33800000 / */
/*     DATA DIVER(1) / 16#34000000 / */
/*     DATA LOG10(1) / 16#3E9A209B / */

/*     MACHINE CONSTANTS FOR THE BURROUGHS 1700 SYSTEM */

/*     DATA RMACH(1) / Z400800000 / */
/*     DATA RMACH(2) / Z5FFFFFFFF / */
/*     DATA RMACH(3) / Z4E9800000 / */
/*     DATA RMACH(4) / Z4EA800000 / */
/*     DATA RMACH(5) / Z500E730E8 / */

/*     MACHINE CONSTANTS FOR THE BURROUGHS 5700/6700/7700 SYSTEMS */

/*     DATA RMACH(1) / O1771000000000000 / */
/*     DATA RMACH(2) / O0777777777777777 / */
/*     DATA RMACH(3) / O1311000000000000 / */
/*     DATA RMACH(4) / O1301000000000000 / */
/*     DATA RMACH(5) / O1157163034761675 / */

/*     MACHINE CONSTANTS FOR THE CDC 170/180 SERIES USING NOS/VE */

/*     DATA RMACH(1) / Z"3001800000000000" / */
/*     DATA RMACH(2) / Z"4FFEFFFFFFFFFFFE" / */
/*     DATA RMACH(3) / Z"3FD2800000000000" / */
/*     DATA RMACH(4) / Z"3FD3800000000000" / */
/*     DATA RMACH(5) / Z"3FFF9A209A84FBCF" / */

/*     MACHINE CONSTANTS FOR THE CDC 6000/7000 SERIES */

/*     DATA RMACH(1) / 00564000000000000000B / */
/*     DATA RMACH(2) / 37767777777777777776B / */
/*     DATA RMACH(3) / 16414000000000000000B / */
/*     DATA RMACH(4) / 16424000000000000000B / */
/*     DATA RMACH(5) / 17164642023241175720B / */

/*     MACHINE CONSTANTS FOR THE CELERITY C1260 */

/*     DATA SMALL(1) / Z'00800000' / */
/*     DATA LARGE(1) / Z'7F7FFFFF' / */
/*     DATA RIGHT(1) / Z'33800000' / */
/*     DATA DIVER(1) / Z'34000000' / */
/*     DATA LOG10(1) / Z'3E9A209B' / */

/*     MACHINE CONSTANTS FOR THE CONVEX */
/*     USING THE -fn COMPILER OPTION */

/*     DATA RMACH(1) / Z'00800000' / */
/*     DATA RMACH(2) / Z'7FFFFFFF' / */
/*     DATA RMACH(3) / Z'34800000' / */
/*     DATA RMACH(4) / Z'35000000' / */
/*     DATA RMACH(5) / Z'3F9A209B' / */

/*     MACHINE CONSTANTS FOR THE CONVEX */
/*     USING THE -fi COMPILER OPTION */

/*     DATA RMACH(1) / Z'00800000' / */
/*     DATA RMACH(2) / Z'7F7FFFFF' / */
/*     DATA RMACH(3) / Z'33800000' / */
/*     DATA RMACH(4) / Z'34000000' / */
/*     DATA RMACH(5) / Z'3E9A209B' / */

/*     MACHINE CONSTANTS FOR THE CONVEX */
/*     USING THE -p8 OR -pd8 COMPILER OPTION */

/*     DATA RMACH(1) / Z'0010000000000000' / */
/*     DATA RMACH(2) / Z'7FFFFFFFFFFFFFFF' / */
/*     DATA RMACH(3) / Z'3CC0000000000000' / */
/*     DATA RMACH(4) / Z'3CD0000000000000' / */
/*     DATA RMACH(5) / Z'3FF34413509F79FF' / */

/*     MACHINE CONSTANTS FOR THE CRAY */

/*     DATA RMACH(1) / 200034000000000000000B / */
/*     DATA RMACH(2) / 577767777777777777776B / */
/*     DATA RMACH(3) / 377224000000000000000B / */
/*     DATA RMACH(4) / 377234000000000000000B / */
/*     DATA RMACH(5) / 377774642023241175720B / */

/*     MACHINE CONSTANTS FOR THE DATA GENERAL ECLIPSE S/200 */
/*     NOTE - IT MAY BE APPROPRIATE TO INCLUDE THE FOLLOWING CARD - */
/*     STATIC RMACH(5) */

/*     DATA SMALL /    20K,       0 / */
/*     DATA LARGE / 77777K, 177777K / */
/*     DATA RIGHT / 35420K,       0 / */
/*     DATA DIVER / 36020K,       0 / */
/*     DATA LOG10 / 40423K,  42023K / */

/*     MACHINE CONSTANTS FOR THE DEC ALPHA */
/*     USING G_FLOAT */

/*     DATA RMACH(1) / '00000080'X / */
/*     DATA RMACH(2) / 'FFFF7FFF'X / */
/*     DATA RMACH(3) / '00003480'X / */
/*     DATA RMACH(4) / '00003500'X / */
/*     DATA RMACH(5) / '209B3F9A'X / */

/*     MACHINE CONSTANTS FOR THE DEC ALPHA */
/*     USING IEEE_FLOAT */

/*     DATA RMACH(1) / '00800000'X / */
/*     DATA RMACH(2) / '7F7FFFFF'X / */
/*     DATA RMACH(3) / '33800000'X / */
/*     DATA RMACH(4) / '34000000'X / */
/*     DATA RMACH(5) / '3E9A209B'X / */

/*     MACHINE CONSTANTS FOR THE DEC RISC */

/*     DATA RMACH(1) / Z'00800000' / */
/*     DATA RMACH(2) / Z'7F7FFFFF' / */
/*     DATA RMACH(3) / Z'33800000' / */
/*     DATA RMACH(4) / Z'34000000' / */
/*     DATA RMACH(5) / Z'3E9A209B' / */

/*     MACHINE CONSTANTS FOR THE DEC VAX */
/*     (EXPRESSED IN INTEGER AND HEXADECIMAL) */
/*     THE HEX FORMAT BELOW MAY NOT BE SUITABLE FOR UNIX SYSTEMS */
/*     THE INTEGER FORMAT SHOULD BE OK FOR UNIX SYSTEMS */

/*     DATA SMALL(1) /       128 / */
/*     DATA LARGE(1) /    -32769 / */
/*     DATA RIGHT(1) /     13440 / */
/*     DATA DIVER(1) /     13568 / */
/*     DATA LOG10(1) / 547045274 / */

/*     DATA SMALL(1) / Z00000080 / */
/*     DATA LARGE(1) / ZFFFF7FFF / */
/*     DATA RIGHT(1) / Z00003480 / */
/*     DATA DIVER(1) / Z00003500 / */
/*     DATA LOG10(1) / Z209B3F9A / */

/*     MACHINE CONSTANTS FOR THE ELXSI 6400 */
/*     (ASSUMING REAL*4 IS THE DEFAULT REAL) */

/*     DATA SMALL(1) / '00800000'X / */
/*     DATA LARGE(1) / '7F7FFFFF'X / */
/*     DATA RIGHT(1) / '33800000'X / */
/*     DATA DIVER(1) / '34000000'X / */
/*     DATA LOG10(1) / '3E9A209B'X / */

/*     MACHINE CONSTANTS FOR THE HARRIS 220 */

/*     DATA SMALL(1), SMALL(2) / '20000000, '00000201 / */
/*     DATA LARGE(1), LARGE(2) / '37777777, '00000177 / */
/*     DATA RIGHT(1), RIGHT(2) / '20000000, '00000352 / */
/*     DATA DIVER(1), DIVER(2) / '20000000, '00000353 / */
/*     DATA LOG10(1), LOG10(2) / '23210115, '00000377 / */

/*     MACHINE CONSTANTS FOR THE HONEYWELL 600/6000 SERIES */

/*     DATA RMACH(1) / O402400000000 / */
/*     DATA RMACH(2) / O376777777777 / */
/*     DATA RMACH(3) / O714400000000 / */
/*     DATA RMACH(4) / O716400000000 / */
/*     DATA RMACH(5) / O776464202324 / */

/*     MACHINE CONSTANTS FOR THE HP 730 */

/*     DATA RMACH(1) / Z'00800000' / */
/*     DATA RMACH(2) / Z'7F7FFFFF' / */
/*     DATA RMACH(3) / Z'33800000' / */
/*     DATA RMACH(4) / Z'34000000' / */
/*     DATA RMACH(5) / Z'3E9A209B' / */

/*     MACHINE CONSTANTS FOR THE HP 2100 */
/*     3 WORD DOUBLE PRECISION WITH FTN4 */

/*     DATA SMALL(1), SMALL(2) / 40000B,       1 / */
/*     DATA LARGE(1), LARGE(2) / 77777B, 177776B / */
/*     DATA RIGHT(1), RIGHT(2) / 40000B,    325B / */
/*     DATA DIVER(1), DIVER(2) / 40000B,    327B / */
/*     DATA LOG10(1), LOG10(2) / 46420B,  46777B / */

/*     MACHINE CONSTANTS FOR THE HP 2100 */
/*     4 WORD DOUBLE PRECISION WITH FTN4 */

/*     DATA SMALL(1), SMALL(2) / 40000B,       1 / */
/*     DATA LARGE(1), LARGE(2) / 77777B, 177776B / */
/*     DATA RIGHT(1), RIGHT(2) / 40000B,    325B / */
/*     DATA DIVER(1), DIVER(2) / 40000B,    327B / */
/*     DATA LOG10(1), LOG10(2) / 46420B,  46777B / */

/*     MACHINE CONSTANTS FOR THE HP 9000 */

/*     DATA SMALL(1) / 00004000000B / */
/*     DATA LARGE(1) / 17677777777B / */
/*     DATA RIGHT(1) / 06340000000B / */
/*     DATA DIVER(1) / 06400000000B / */
/*     DATA LOG10(1) / 07646420233B / */

/*     MACHINE CONSTANTS FOR THE IBM 360/370 SERIES, */
/*     THE XEROX SIGMA 5/7/9, THE SEL SYSTEMS 85/86  AND */
/*     THE PERKIN ELMER (INTERDATA) 7/32. */

/*     DATA RMACH(1) / Z00100000 / */
/*     DATA RMACH(2) / Z7FFFFFFF / */
/*     DATA RMACH(3) / Z3B100000 / */
/*     DATA RMACH(4) / Z3C100000 / */
/*     DATA RMACH(5) / Z41134413 / */

/*     MACHINE CONSTANTS FOR THE IBM PC */

/*     DATA SMALL(1) / 1.18E-38      / */
/*     DATA LARGE(1) / 3.40E+38      / */
/*     DATA RIGHT(1) / 0.595E-07     / */
/*     DATA DIVER(1) / 1.19E-07      / */
/*     DATA LOG10(1) / 0.30102999566 / */

/*     MACHINE CONSTANTS FOR THE IBM RS 6000 */

/*     DATA RMACH(1) / Z'00800000' / */
/*     DATA RMACH(2) / Z'7F7FFFFF' / */
/*     DATA RMACH(3) / Z'33800000' / */
/*     DATA RMACH(4) / Z'34000000' / */
/*     DATA RMACH(5) / Z'3E9A209B' / */

/*     MACHINE CONSTANTS FOR THE INTEL i860 */

/*     DATA RMACH(1) / Z'00800000' / */
/*     DATA RMACH(2) / Z'7F7FFFFF' / */
/*     DATA RMACH(3) / Z'33800000' / */
/*     DATA RMACH(4) / Z'34000000' / */
/*     DATA RMACH(5) / Z'3E9A209B' / */

/*     MACHINE CONSTANTS FOR THE PDP-10 (KA OR KI PROCESSOR) */

/*     DATA RMACH(1) / "000400000000 / */
/*     DATA RMACH(2) / "377777777777 / */
/*     DATA RMACH(3) / "146400000000 / */
/*     DATA RMACH(4) / "147400000000 / */
/*     DATA RMACH(5) / "177464202324 / */

/*     MACHINE CONSTANTS FOR PDP-11 FORTRAN SUPPORTING */
/*     32-BIT INTEGERS (EXPRESSED IN INTEGER AND OCTAL). */

/*     DATA SMALL(1) /    8388608 / */
/*     DATA LARGE(1) / 2147483647 / */
/*     DATA RIGHT(1) /  880803840 / */
/*     DATA DIVER(1) /  889192448 / */
/*     DATA LOG10(1) / 1067065499 / */

/*     DATA RMACH(1) / O00040000000 / */
/*     DATA RMACH(2) / O17777777777 / */
/*     DATA RMACH(3) / O06440000000 / */
/*     DATA RMACH(4) / O06500000000 / */
/*     DATA RMACH(5) / O07746420233 / */

/*     MACHINE CONSTANTS FOR PDP-11 FORTRAN SUPPORTING */
/*     16-BIT INTEGERS  (EXPRESSED IN INTEGER AND OCTAL). */

/*     DATA SMALL(1), SMALL(2) /   128,     0 / */
/*     DATA LARGE(1), LARGE(2) / 32767,    -1 / */
/*     DATA RIGHT(1), RIGHT(2) / 13440,     0 / */
/*     DATA DIVER(1), DIVER(2) / 13568,     0 / */
/*     DATA LOG10(1), LOG10(2) / 16282,  8347 / */

/*     DATA SMALL(1), SMALL(2) / O000200, O000000 / */
/*     DATA LARGE(1), LARGE(2) / O077777, O177777 / */
/*     DATA RIGHT(1), RIGHT(2) / O032200, O000000 / */
/*     DATA DIVER(1), DIVER(2) / O032400, O000000 / */
/*     DATA LOG10(1), LOG10(2) / O037632, O020233 / */

/*     MACHINE CONSTANTS FOR THE SILICON GRAPHICS */

/*     DATA RMACH(1) / Z'00800000' / */
/*     DATA RMACH(2) / Z'7F7FFFFF' / */
/*     DATA RMACH(3) / Z'33800000' / */
/*     DATA RMACH(4) / Z'34000000' / */
/*     DATA RMACH(5) / Z'3E9A209B' / */

/*     MACHINE CONSTANTS FOR THE SUN */

/*     DATA RMACH(1) / Z'00800000' / */
/*     DATA RMACH(2) / Z'7F7FFFFF' / */
/*     DATA RMACH(3) / Z'33800000' / */
/*     DATA RMACH(4) / Z'34000000' / */
/*     DATA RMACH(5) / Z'3E9A209B' / */

/*     MACHINE CONSTANTS FOR THE SUN */
/*     USING THE -r8 COMPILER OPTION */

/*     DATA RMACH(1) / Z'0010000000000000' / */
/*     DATA RMACH(2) / Z'7FEFFFFFFFFFFFFF' / */
/*     DATA RMACH(3) / Z'3CA0000000000000' / */
/*     DATA RMACH(4) / Z'3CB0000000000000' / */
/*     DATA RMACH(5) / Z'3FD34413509F79FF' / */

/*     MACHINE CONSTANTS FOR THE UNIVAC 1100 SERIES */

/*     DATA RMACH(1) / O000400000000 / */
/*     DATA RMACH(2) / O377777777777 / */
/*     DATA RMACH(3) / O146400000000 / */
/*     DATA RMACH(4) / O147400000000 / */
/*     DATA RMACH(5) / O177464202324 / */

/*     MACHINE CONSTANTS FOR THE Z80 MICROPROCESSOR */

/*     DATA SMALL(1), SMALL(2) /     0,    256/ */
/*     DATA LARGE(1), LARGE(2) /    -1,   -129/ */
/*     DATA RIGHT(1), RIGHT(2) /     0,  26880/ */
/*     DATA DIVER(1), DIVER(2) /     0,  27136/ */
/*     DATA LOG10(1), LOG10(2) /  8347,  32538/ */

/* ***FIRST EXECUTABLE STATEMENT  R1MACH */
    if (*i__ < 1 || *i__ > 5) {
	xermsg_("SLATEC", "R1MACH", "I OUT OF BOUNDS", &c__1, &c__2, (ftnlen)
		6, (ftnlen)6, (ftnlen)15);
    }

    ret_val = rmach[*i__ - 1];
    return ret_val;

} /* r1mach_ */

#undef right
#undef diver
#undef small
#undef rmach
#undef large
#undef log10


