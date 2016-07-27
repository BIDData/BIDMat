/* i1mach.f -- translated by f2c (version 20100827).
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
__device__ __constant__ integer equiv_0[16] = {
  5,
  6,
  6,
  6,
  32,
  4,
  2,
  31,
  2147483647,
  2,
  24,
  -125,
  128,
  53,
  -1021,
  1024};

/* DECK I1MACH */
__device__ integer i1mach_(integer *i__)
{
    /* Format strings */
    const char fmt_9000[] = "(\0021ERROR    1 IN I1MACH - I OUT OF BOUND"
	    "S\002)";

    /* System generated locals */
    integer ret_val;

    /* Builtin functions */
    integer s_wsfe(cilist *), e_wsfe(void);
    /* Subroutine */ int s_stop(char *, ftnlen);

    /* Local variables */
#define imach (equiv_0)
#define output (equiv_0 + 3)

    /* Fortran I/O blocks */
    //    static cilist io___3 = { 0, 0, 0, fmt_9000, 0 };


/* ***BEGIN PROLOGUE  I1MACH */
/* ***PURPOSE  Return integer machine dependent constants. */
/* ***LIBRARY   SLATEC */
/* ***CATEGORY  R1 */
/* ***TYPE      INTEGER (I1MACH-I) */
/* ***KEYWORDS  MACHINE CONSTANTS */
/* ***AUTHOR  Fox, P. A., (Bell Labs) */
/*           Hall, A. D., (Bell Labs) */
/*           Schryer, N. L., (Bell Labs) */
/* ***DESCRIPTION */

/*   I1MACH can be used to obtain machine-dependent parameters for the */
/*   local machine environment.  It is a function subprogram with one */
/*   (input) argument and can be referenced as follows: */

/*        K = I1MACH(I) */

/*   where I=1,...,16.  The (output) value of K above is determined by */
/*   the (input) value of I.  The results for various values of I are */
/*   discussed below. */

/*   I/O unit numbers: */
/*     I1MACH( 1) = the standard input unit. */
/*     I1MACH( 2) = the standard output unit. */
/*     I1MACH( 3) = the standard punch unit. */
/*     I1MACH( 4) = the standard error message unit. */

/*   Words: */
/*     I1MACH( 5) = the number of bits per integer storage unit. */
/*     I1MACH( 6) = the number of characters per integer storage unit. */

/*   Integers: */
/*     assume integers are represented in the S-digit, base-A form */

/*                sign ( X(S-1)*A**(S-1) + ... + X(1)*A + X(0) ) */

/*                where 0 .LE. X(I) .LT. A for I=0,...,S-1. */
/*     I1MACH( 7) = A, the base. */
/*     I1MACH( 8) = S, the number of base-A digits. */
/*     I1MACH( 9) = A**S - 1, the largest magnitude. */

/*   Floating-Point Numbers: */
/*     Assume floating-point numbers are represented in the T-digit, */
/*     base-B form */
/*                sign (B**E)*( (X(1)/B) + ... + (X(T)/B**T) ) */

/*                where 0 .LE. X(I) .LT. B for I=1,...,T, */
/*                0 .LT. X(1), and EMIN .LE. E .LE. EMAX. */
/*     I1MACH(10) = B, the base. */

/*   Single-Precision: */
/*     I1MACH(11) = T, the number of base-B digits. */
/*     I1MACH(12) = EMIN, the smallest exponent E. */
/*     I1MACH(13) = EMAX, the largest exponent E. */

/*   Double-Precision: */
/*     I1MACH(14) = T, the number of base-B digits. */
/*     I1MACH(15) = EMIN, the smallest exponent E. */
/*     I1MACH(16) = EMAX, the largest exponent E. */

/*   To alter this function for a particular environment, the desired */
/*   set of DATA statements should be activated by removing the C from */
/*   column 1.  Also, the values of I1MACH(1) - I1MACH(4) should be */
/*   checked for consistency with the local operating system. */

/* ***REFERENCES  P. A. Fox, A. D. Hall and N. L. Schryer, Framework for */
/*                 a portable library, ACM Transactions on Mathematical */
/*                 Software 4, 2 (June 1978), pp. 177-188. */
/* ***ROUTINES CALLED  (NONE) */
/* ***REVISION HISTORY  (YYMMDD) */
/*   750101  DATE WRITTEN */
/*   891012  Added VAX G-floating constants.  (WRB) */
/*   891012  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900618  Added DEC RISC constants.  (WRB) */
/*   900723  Added IBM RS 6000 constants.  (WRB) */
/*   901009  Correct I1MACH(7) for IBM Mainframes. Should be 2 not 16. */
/*           (RWC) */
/*   910710  Added HP 730 constants.  (SMR) */
/*   911114  Added Convex IEEE constants.  (WRB) */
/*   920121  Added SUN -r8 compiler option constants.  (WRB) */
/*   920229  Added Touchstone Delta i860 constants.  (WRB) */
/*   920501  Reformatted the REFERENCES section.  (WRB) */
/*   920625  Added Convex -p8 and -pd8 compiler option constants. */
/*           (BKS, WRB) */
/*   930201  Added DEC Alpha and SGI constants.  (RWC and WRB) */
/*   930618  Corrected I1MACH(5) for Convex -p8 and -pd8 compiler */
/*           options.  (DWL, RWC and WRB). */
/* ***END PROLOGUE  I1MACH */


/*     MACHINE CONSTANTS FOR THE AMIGA */
/*     ABSOFT COMPILER */

/*     DATA IMACH( 1) /          5 / */
/*     DATA IMACH( 2) /          6 / */
/*     DATA IMACH( 3) /          5 / */
/*     DATA IMACH( 4) /          6 / */
/*     DATA IMACH( 5) /         32 / */
/*     DATA IMACH( 6) /          4 / */
/*     DATA IMACH( 7) /          2 / */
/*     DATA IMACH( 8) /         31 / */
/*     DATA IMACH( 9) / 2147483647 / */
/*     DATA IMACH(10) /          2 / */
/*     DATA IMACH(11) /         24 / */
/*     DATA IMACH(12) /       -126 / */
/*     DATA IMACH(13) /        127 / */
/*     DATA IMACH(14) /         53 / */
/*     DATA IMACH(15) /      -1022 / */
/*     DATA IMACH(16) /       1023 / */

/*     MACHINE CONSTANTS FOR THE APOLLO */

/*     DATA IMACH( 1) /          5 / */
/*     DATA IMACH( 2) /          6 / */
/*     DATA IMACH( 3) /          6 / */
/*     DATA IMACH( 4) /          6 / */
/*     DATA IMACH( 5) /         32 / */
/*     DATA IMACH( 6) /          4 / */
/*     DATA IMACH( 7) /          2 / */
/*     DATA IMACH( 8) /         31 / */
/*     DATA IMACH( 9) / 2147483647 / */
/*     DATA IMACH(10) /          2 / */
/*     DATA IMACH(11) /         24 / */
/*     DATA IMACH(12) /       -125 / */
/*     DATA IMACH(13) /        129 / */
/*     DATA IMACH(14) /         53 / */
/*     DATA IMACH(15) /      -1021 / */
/*     DATA IMACH(16) /       1025 / */

/*     MACHINE CONSTANTS FOR THE BURROUGHS 1700 SYSTEM */

/*     DATA IMACH( 1) /          7 / */
/*     DATA IMACH( 2) /          2 / */
/*     DATA IMACH( 3) /          2 / */
/*     DATA IMACH( 4) /          2 / */
/*     DATA IMACH( 5) /         36 / */
/*     DATA IMACH( 6) /          4 / */
/*     DATA IMACH( 7) /          2 / */
/*     DATA IMACH( 8) /         33 / */
/*     DATA IMACH( 9) / Z1FFFFFFFF / */
/*     DATA IMACH(10) /          2 / */
/*     DATA IMACH(11) /         24 / */
/*     DATA IMACH(12) /       -256 / */
/*     DATA IMACH(13) /        255 / */
/*     DATA IMACH(14) /         60 / */
/*     DATA IMACH(15) /       -256 / */
/*     DATA IMACH(16) /        255 / */

/*     MACHINE CONSTANTS FOR THE BURROUGHS 5700 SYSTEM */

/*     DATA IMACH( 1) /          5 / */
/*     DATA IMACH( 2) /          6 / */
/*     DATA IMACH( 3) /          7 / */
/*     DATA IMACH( 4) /          6 / */
/*     DATA IMACH( 5) /         48 / */
/*     DATA IMACH( 6) /          6 / */
/*     DATA IMACH( 7) /          2 / */
/*     DATA IMACH( 8) /         39 / */
/*     DATA IMACH( 9) / O0007777777777777 / */
/*     DATA IMACH(10) /          8 / */
/*     DATA IMACH(11) /         13 / */
/*     DATA IMACH(12) /        -50 / */
/*     DATA IMACH(13) /         76 / */
/*     DATA IMACH(14) /         26 / */
/*     DATA IMACH(15) /        -50 / */
/*     DATA IMACH(16) /         76 / */

/*     MACHINE CONSTANTS FOR THE BURROUGHS 6700/7700 SYSTEMS */

/*     DATA IMACH( 1) /          5 / */
/*     DATA IMACH( 2) /          6 / */
/*     DATA IMACH( 3) /          7 / */
/*     DATA IMACH( 4) /          6 / */
/*     DATA IMACH( 5) /         48 / */
/*     DATA IMACH( 6) /          6 / */
/*     DATA IMACH( 7) /          2 / */
/*     DATA IMACH( 8) /         39 / */
/*     DATA IMACH( 9) / O0007777777777777 / */
/*     DATA IMACH(10) /          8 / */
/*     DATA IMACH(11) /         13 / */
/*     DATA IMACH(12) /        -50 / */
/*     DATA IMACH(13) /         76 / */
/*     DATA IMACH(14) /         26 / */
/*     DATA IMACH(15) /     -32754 / */
/*     DATA IMACH(16) /      32780 / */

/*     MACHINE CONSTANTS FOR THE CDC 170/180 SERIES USING NOS/VE */

/*     DATA IMACH( 1) /          5 / */
/*     DATA IMACH( 2) /          6 / */
/*     DATA IMACH( 3) /          7 / */
/*     DATA IMACH( 4) /          6 / */
/*     DATA IMACH( 5) /         64 / */
/*     DATA IMACH( 6) /          8 / */
/*     DATA IMACH( 7) /          2 / */
/*     DATA IMACH( 8) /         63 / */
/*     DATA IMACH( 9) / 9223372036854775807 / */
/*     DATA IMACH(10) /          2 / */
/*     DATA IMACH(11) /         47 / */
/*     DATA IMACH(12) /      -4095 / */
/*     DATA IMACH(13) /       4094 / */
/*     DATA IMACH(14) /         94 / */
/*     DATA IMACH(15) /      -4095 / */
/*     DATA IMACH(16) /       4094 / */

/*     MACHINE CONSTANTS FOR THE CDC 6000/7000 SERIES */

/*     DATA IMACH( 1) /          5 / */
/*     DATA IMACH( 2) /          6 / */
/*     DATA IMACH( 3) /          7 / */
/*     DATA IMACH( 4) /    6LOUTPUT/ */
/*     DATA IMACH( 5) /         60 / */
/*     DATA IMACH( 6) /         10 / */
/*     DATA IMACH( 7) /          2 / */
/*     DATA IMACH( 8) /         48 / */
/*     DATA IMACH( 9) / 00007777777777777777B / */
/*     DATA IMACH(10) /          2 / */
/*     DATA IMACH(11) /         47 / */
/*     DATA IMACH(12) /       -929 / */
/*     DATA IMACH(13) /       1070 / */
/*     DATA IMACH(14) /         94 / */
/*     DATA IMACH(15) /       -929 / */
/*     DATA IMACH(16) /       1069 / */

/*     MACHINE CONSTANTS FOR THE CELERITY C1260 */

/*     DATA IMACH( 1) /          5 / */
/*     DATA IMACH( 2) /          6 / */
/*     DATA IMACH( 3) /          6 / */
/*     DATA IMACH( 4) /          0 / */
/*     DATA IMACH( 5) /         32 / */
/*     DATA IMACH( 6) /          4 / */
/*     DATA IMACH( 7) /          2 / */
/*     DATA IMACH( 8) /         31 / */
/*     DATA IMACH( 9) / Z'7FFFFFFF' / */
/*     DATA IMACH(10) /          2 / */
/*     DATA IMACH(11) /         24 / */
/*     DATA IMACH(12) /       -126 / */
/*     DATA IMACH(13) /        127 / */
/*     DATA IMACH(14) /         53 / */
/*     DATA IMACH(15) /      -1022 / */
/*     DATA IMACH(16) /       1023 / */

/*     MACHINE CONSTANTS FOR THE CONVEX */
/*     USING THE -fn COMPILER OPTION */

/*     DATA IMACH( 1) /          5 / */
/*     DATA IMACH( 2) /          6 / */
/*     DATA IMACH( 3) /          7 / */
/*     DATA IMACH( 4) /          6 / */
/*     DATA IMACH( 5) /         32 / */
/*     DATA IMACH( 6) /          4 / */
/*     DATA IMACH( 7) /          2 / */
/*     DATA IMACH( 8) /         31 / */
/*     DATA IMACH( 9) / 2147483647 / */
/*     DATA IMACH(10) /          2 / */
/*     DATA IMACH(11) /         24 / */
/*     DATA IMACH(12) /       -127 / */
/*     DATA IMACH(13) /        127 / */
/*     DATA IMACH(14) /         53 / */
/*     DATA IMACH(15) /      -1023 / */
/*     DATA IMACH(16) /       1023 / */

/*     MACHINE CONSTANTS FOR THE CONVEX */
/*     USING THE -fi COMPILER OPTION */

/*     DATA IMACH( 1) /          5 / */
/*     DATA IMACH( 2) /          6 / */
/*     DATA IMACH( 3) /          7 / */
/*     DATA IMACH( 4) /          6 / */
/*     DATA IMACH( 5) /         32 / */
/*     DATA IMACH( 6) /          4 / */
/*     DATA IMACH( 7) /          2 / */
/*     DATA IMACH( 8) /         31 / */
/*     DATA IMACH( 9) / 2147483647 / */
/*     DATA IMACH(10) /          2 / */
/*     DATA IMACH(11) /         24 / */
/*     DATA IMACH(12) /       -125 / */
/*     DATA IMACH(13) /        128 / */
/*     DATA IMACH(14) /         53 / */
/*     DATA IMACH(15) /      -1021 / */
/*     DATA IMACH(16) /       1024 / */

/*     MACHINE CONSTANTS FOR THE CONVEX */
/*     USING THE -p8 COMPILER OPTION */

/*     DATA IMACH( 1) /          5 / */
/*     DATA IMACH( 2) /          6 / */
/*     DATA IMACH( 3) /          7 / */
/*     DATA IMACH( 4) /          6 / */
/*     DATA IMACH( 5) /         64 / */
/*     DATA IMACH( 6) /          4 / */
/*     DATA IMACH( 7) /          2 / */
/*     DATA IMACH( 8) /         63 / */
/*     DATA IMACH( 9) / 9223372036854775807 / */
/*     DATA IMACH(10) /          2 / */
/*     DATA IMACH(11) /         53 / */
/*     DATA IMACH(12) /      -1023 / */
/*     DATA IMACH(13) /       1023 / */
/*     DATA IMACH(14) /        113 / */
/*     DATA IMACH(15) /     -16383 / */
/*     DATA IMACH(16) /      16383 / */

/*     MACHINE CONSTANTS FOR THE CONVEX */
/*     USING THE -pd8 COMPILER OPTION */

/*     DATA IMACH( 1) /          5 / */
/*     DATA IMACH( 2) /          6 / */
/*     DATA IMACH( 3) /          7 / */
/*     DATA IMACH( 4) /          6 / */
/*     DATA IMACH( 5) /         64 / */
/*     DATA IMACH( 6) /          4 / */
/*     DATA IMACH( 7) /          2 / */
/*     DATA IMACH( 8) /         63 / */
/*     DATA IMACH( 9) / 9223372036854775807 / */
/*     DATA IMACH(10) /          2 / */
/*     DATA IMACH(11) /         53 / */
/*     DATA IMACH(12) /      -1023 / */
/*     DATA IMACH(13) /       1023 / */
/*     DATA IMACH(14) /         53 / */
/*     DATA IMACH(15) /      -1023 / */
/*     DATA IMACH(16) /       1023 / */

/*     MACHINE CONSTANTS FOR THE CRAY */
/*     USING THE 46 BIT INTEGER COMPILER OPTION */

/*     DATA IMACH( 1) /        100 / */
/*     DATA IMACH( 2) /        101 / */
/*     DATA IMACH( 3) /        102 / */
/*     DATA IMACH( 4) /        101 / */
/*     DATA IMACH( 5) /         64 / */
/*     DATA IMACH( 6) /          8 / */
/*     DATA IMACH( 7) /          2 / */
/*     DATA IMACH( 8) /         46 / */
/*     DATA IMACH( 9) / 1777777777777777B / */
/*     DATA IMACH(10) /          2 / */
/*     DATA IMACH(11) /         47 / */
/*     DATA IMACH(12) /      -8189 / */
/*     DATA IMACH(13) /       8190 / */
/*     DATA IMACH(14) /         94 / */
/*     DATA IMACH(15) /      -8099 / */
/*     DATA IMACH(16) /       8190 / */

/*     MACHINE CONSTANTS FOR THE CRAY */
/*     USING THE 64 BIT INTEGER COMPILER OPTION */

/*     DATA IMACH( 1) /        100 / */
/*     DATA IMACH( 2) /        101 / */
/*     DATA IMACH( 3) /        102 / */
/*     DATA IMACH( 4) /        101 / */
/*     DATA IMACH( 5) /         64 / */
/*     DATA IMACH( 6) /          8 / */
/*     DATA IMACH( 7) /          2 / */
/*     DATA IMACH( 8) /         63 / */
/*     DATA IMACH( 9) / 777777777777777777777B / */
/*     DATA IMACH(10) /          2 / */
/*     DATA IMACH(11) /         47 / */
/*     DATA IMACH(12) /      -8189 / */
/*     DATA IMACH(13) /       8190 / */
/*     DATA IMACH(14) /         94 / */
/*     DATA IMACH(15) /      -8099 / */
/*     DATA IMACH(16) /       8190 / */

/*     MACHINE CONSTANTS FOR THE DATA GENERAL ECLIPSE S/200 */

/*     DATA IMACH( 1) /         11 / */
/*     DATA IMACH( 2) /         12 / */
/*     DATA IMACH( 3) /          8 / */
/*     DATA IMACH( 4) /         10 / */
/*     DATA IMACH( 5) /         16 / */
/*     DATA IMACH( 6) /          2 / */
/*     DATA IMACH( 7) /          2 / */
/*     DATA IMACH( 8) /         15 / */
/*     DATA IMACH( 9) /      32767 / */
/*     DATA IMACH(10) /         16 / */
/*     DATA IMACH(11) /          6 / */
/*     DATA IMACH(12) /        -64 / */
/*     DATA IMACH(13) /         63 / */
/*     DATA IMACH(14) /         14 / */
/*     DATA IMACH(15) /        -64 / */
/*     DATA IMACH(16) /         63 / */

/*     MACHINE CONSTANTS FOR THE DEC ALPHA */
/*     USING G_FLOAT */

/*     DATA IMACH( 1) /          5 / */
/*     DATA IMACH( 2) /          6 / */
/*     DATA IMACH( 3) /          5 / */
/*     DATA IMACH( 4) /          6 / */
/*     DATA IMACH( 5) /         32 / */
/*     DATA IMACH( 6) /          4 / */
/*     DATA IMACH( 7) /          2 / */
/*     DATA IMACH( 8) /         31 / */
/*     DATA IMACH( 9) / 2147483647 / */
/*     DATA IMACH(10) /          2 / */
/*     DATA IMACH(11) /         24 / */
/*     DATA IMACH(12) /       -127 / */
/*     DATA IMACH(13) /        127 / */
/*     DATA IMACH(14) /         53 / */
/*     DATA IMACH(15) /      -1023 / */
/*     DATA IMACH(16) /       1023 / */

/*     MACHINE CONSTANTS FOR THE DEC ALPHA */
/*     USING IEEE_FLOAT */

/*     DATA IMACH( 1) /          5 / */
/*     DATA IMACH( 2) /          6 / */
/*     DATA IMACH( 3) /          6 / */
/*     DATA IMACH( 4) /          6 / */
/*     DATA IMACH( 5) /         32 / */
/*     DATA IMACH( 6) /          4 / */
/*     DATA IMACH( 7) /          2 / */
/*     DATA IMACH( 8) /         31 / */
/*     DATA IMACH( 9) / 2147483647 / */
/*     DATA IMACH(10) /          2 / */
/*     DATA IMACH(11) /         24 / */
/*     DATA IMACH(12) /       -125 / */
/*     DATA IMACH(13) /        128 / */
/*     DATA IMACH(14) /         53 / */
/*     DATA IMACH(15) /      -1021 / */
/*     DATA IMACH(16) /       1024 / */

/*     MACHINE CONSTANTS FOR THE DEC RISC */

/*     DATA IMACH( 1) /          5 / */
/*     DATA IMACH( 2) /          6 / */
/*     DATA IMACH( 3) /          6 / */
/*     DATA IMACH( 4) /          6 / */
/*     DATA IMACH( 5) /         32 / */
/*     DATA IMACH( 6) /          4 / */
/*     DATA IMACH( 7) /          2 / */
/*     DATA IMACH( 8) /         31 / */
/*     DATA IMACH( 9) / 2147483647 / */
/*     DATA IMACH(10) /          2 / */
/*     DATA IMACH(11) /         24 / */
/*     DATA IMACH(12) /       -125 / */
/*     DATA IMACH(13) /        128 / */
/*     DATA IMACH(14) /         53 / */
/*     DATA IMACH(15) /      -1021 / */
/*     DATA IMACH(16) /       1024 / */

/*     MACHINE CONSTANTS FOR THE DEC VAX */
/*     USING D_FLOATING */

/*     DATA IMACH( 1) /          5 / */
/*     DATA IMACH( 2) /          6 / */
/*     DATA IMACH( 3) /          5 / */
/*     DATA IMACH( 4) /          6 / */
/*     DATA IMACH( 5) /         32 / */
/*     DATA IMACH( 6) /          4 / */
/*     DATA IMACH( 7) /          2 / */
/*     DATA IMACH( 8) /         31 / */
/*     DATA IMACH( 9) / 2147483647 / */
/*     DATA IMACH(10) /          2 / */
/*     DATA IMACH(11) /         24 / */
/*     DATA IMACH(12) /       -127 / */
/*     DATA IMACH(13) /        127 / */
/*     DATA IMACH(14) /         56 / */
/*     DATA IMACH(15) /       -127 / */
/*     DATA IMACH(16) /        127 / */

/*     MACHINE CONSTANTS FOR THE DEC VAX */
/*     USING G_FLOATING */

/*     DATA IMACH( 1) /          5 / */
/*     DATA IMACH( 2) /          6 / */
/*     DATA IMACH( 3) /          5 / */
/*     DATA IMACH( 4) /          6 / */
/*     DATA IMACH( 5) /         32 / */
/*     DATA IMACH( 6) /          4 / */
/*     DATA IMACH( 7) /          2 / */
/*     DATA IMACH( 8) /         31 / */
/*     DATA IMACH( 9) / 2147483647 / */
/*     DATA IMACH(10) /          2 / */
/*     DATA IMACH(11) /         24 / */
/*     DATA IMACH(12) /       -127 / */
/*     DATA IMACH(13) /        127 / */
/*     DATA IMACH(14) /         53 / */
/*     DATA IMACH(15) /      -1023 / */
/*     DATA IMACH(16) /       1023 / */

/*     MACHINE CONSTANTS FOR THE ELXSI 6400 */

/*     DATA IMACH( 1) /          5 / */
/*     DATA IMACH( 2) /          6 / */
/*     DATA IMACH( 3) /          6 / */
/*     DATA IMACH( 4) /          6 / */
/*     DATA IMACH( 5) /         32 / */
/*     DATA IMACH( 6) /          4 / */
/*     DATA IMACH( 7) /          2 / */
/*     DATA IMACH( 8) /         32 / */
/*     DATA IMACH( 9) / 2147483647 / */
/*     DATA IMACH(10) /          2 / */
/*     DATA IMACH(11) /         24 / */
/*     DATA IMACH(12) /       -126 / */
/*     DATA IMACH(13) /        127 / */
/*     DATA IMACH(14) /         53 / */
/*     DATA IMACH(15) /      -1022 / */
/*     DATA IMACH(16) /       1023 / */

/*     MACHINE CONSTANTS FOR THE HARRIS 220 */

/*     DATA IMACH( 1) /          5 / */
/*     DATA IMACH( 2) /          6 / */
/*     DATA IMACH( 3) /          0 / */
/*     DATA IMACH( 4) /          6 / */
/*     DATA IMACH( 5) /         24 / */
/*     DATA IMACH( 6) /          3 / */
/*     DATA IMACH( 7) /          2 / */
/*     DATA IMACH( 8) /         23 / */
/*     DATA IMACH( 9) /    8388607 / */
/*     DATA IMACH(10) /          2 / */
/*     DATA IMACH(11) /         23 / */
/*     DATA IMACH(12) /       -127 / */
/*     DATA IMACH(13) /        127 / */
/*     DATA IMACH(14) /         38 / */
/*     DATA IMACH(15) /       -127 / */
/*     DATA IMACH(16) /        127 / */

/*     MACHINE CONSTANTS FOR THE HONEYWELL 600/6000 SERIES */

/*     DATA IMACH( 1) /          5 / */
/*     DATA IMACH( 2) /          6 / */
/*     DATA IMACH( 3) /         43 / */
/*     DATA IMACH( 4) /          6 / */
/*     DATA IMACH( 5) /         36 / */
/*     DATA IMACH( 6) /          6 / */
/*     DATA IMACH( 7) /          2 / */
/*     DATA IMACH( 8) /         35 / */
/*     DATA IMACH( 9) / O377777777777 / */
/*     DATA IMACH(10) /          2 / */
/*     DATA IMACH(11) /         27 / */
/*     DATA IMACH(12) /       -127 / */
/*     DATA IMACH(13) /        127 / */
/*     DATA IMACH(14) /         63 / */
/*     DATA IMACH(15) /       -127 / */
/*     DATA IMACH(16) /        127 / */

/*     MACHINE CONSTANTS FOR THE HP 730 */

/*     DATA IMACH( 1) /          5 / */
/*     DATA IMACH( 2) /          6 / */
/*     DATA IMACH( 3) /          6 / */
/*     DATA IMACH( 4) /          6 / */
/*     DATA IMACH( 5) /         32 / */
/*     DATA IMACH( 6) /          4 / */
/*     DATA IMACH( 7) /          2 / */
/*     DATA IMACH( 8) /         31 / */
/*     DATA IMACH( 9) / 2147483647 / */
/*     DATA IMACH(10) /          2 / */
/*     DATA IMACH(11) /         24 / */
/*     DATA IMACH(12) /       -125 / */
/*     DATA IMACH(13) /        128 / */
/*     DATA IMACH(14) /         53 / */
/*     DATA IMACH(15) /      -1021 / */
/*     DATA IMACH(16) /       1024 / */

/*     MACHINE CONSTANTS FOR THE HP 2100 */
/*     3 WORD DOUBLE PRECISION OPTION WITH FTN4 */

/*     DATA IMACH( 1) /          5 / */
/*     DATA IMACH( 2) /          6 / */
/*     DATA IMACH( 3) /          4 / */
/*     DATA IMACH( 4) /          1 / */
/*     DATA IMACH( 5) /         16 / */
/*     DATA IMACH( 6) /          2 / */
/*     DATA IMACH( 7) /          2 / */
/*     DATA IMACH( 8) /         15 / */
/*     DATA IMACH( 9) /      32767 / */
/*     DATA IMACH(10) /          2 / */
/*     DATA IMACH(11) /         23 / */
/*     DATA IMACH(12) /       -128 / */
/*     DATA IMACH(13) /        127 / */
/*     DATA IMACH(14) /         39 / */
/*     DATA IMACH(15) /       -128 / */
/*     DATA IMACH(16) /        127 / */

/*     MACHINE CONSTANTS FOR THE HP 2100 */
/*     4 WORD DOUBLE PRECISION OPTION WITH FTN4 */

/*     DATA IMACH( 1) /          5 / */
/*     DATA IMACH( 2) /          6 / */
/*     DATA IMACH( 3) /          4 / */
/*     DATA IMACH( 4) /          1 / */
/*     DATA IMACH( 5) /         16 / */
/*     DATA IMACH( 6) /          2 / */
/*     DATA IMACH( 7) /          2 / */
/*     DATA IMACH( 8) /         15 / */
/*     DATA IMACH( 9) /      32767 / */
/*     DATA IMACH(10) /          2 / */
/*     DATA IMACH(11) /         23 / */
/*     DATA IMACH(12) /       -128 / */
/*     DATA IMACH(13) /        127 / */
/*     DATA IMACH(14) /         55 / */
/*     DATA IMACH(15) /       -128 / */
/*     DATA IMACH(16) /        127 / */

/*     MACHINE CONSTANTS FOR THE HP 9000 */

/*     DATA IMACH( 1) /          5 / */
/*     DATA IMACH( 2) /          6 / */
/*     DATA IMACH( 3) /          6 / */
/*     DATA IMACH( 4) /          7 / */
/*     DATA IMACH( 5) /         32 / */
/*     DATA IMACH( 6) /          4 / */
/*     DATA IMACH( 7) /          2 / */
/*     DATA IMACH( 8) /         32 / */
/*     DATA IMACH( 9) / 2147483647 / */
/*     DATA IMACH(10) /          2 / */
/*     DATA IMACH(11) /         24 / */
/*     DATA IMACH(12) /       -126 / */
/*     DATA IMACH(13) /        127 / */
/*     DATA IMACH(14) /         53 / */
/*     DATA IMACH(15) /      -1015 / */
/*     DATA IMACH(16) /       1017 / */

/*     MACHINE CONSTANTS FOR THE IBM 360/370 SERIES, */
/*     THE XEROX SIGMA 5/7/9, THE SEL SYSTEMS 85/86, AND */
/*     THE PERKIN ELMER (INTERDATA) 7/32. */

/*     DATA IMACH( 1) /          5 / */
/*     DATA IMACH( 2) /          6 / */
/*     DATA IMACH( 3) /          7 / */
/*     DATA IMACH( 4) /          6 / */
/*     DATA IMACH( 5) /         32 / */
/*     DATA IMACH( 6) /          4 / */
/*     DATA IMACH( 7) /          2 / */
/*     DATA IMACH( 8) /         31 / */
/*     DATA IMACH( 9) /  Z7FFFFFFF / */
/*     DATA IMACH(10) /         16 / */
/*     DATA IMACH(11) /          6 / */
/*     DATA IMACH(12) /        -64 / */
/*     DATA IMACH(13) /         63 / */
/*     DATA IMACH(14) /         14 / */
/*     DATA IMACH(15) /        -64 / */
/*     DATA IMACH(16) /         63 / */

/*     MACHINE CONSTANTS FOR THE IBM PC */

/*     DATA IMACH( 1) /          5 / */
/*     DATA IMACH( 2) /          6 / */
/*     DATA IMACH( 3) /          0 / */
/*     DATA IMACH( 4) /          0 / */
/*     DATA IMACH( 5) /         32 / */
/*     DATA IMACH( 6) /          4 / */
/*     DATA IMACH( 7) /          2 / */
/*     DATA IMACH( 8) /         31 / */
/*     DATA IMACH( 9) / 2147483647 / */
/*     DATA IMACH(10) /          2 / */
/*     DATA IMACH(11) /         24 / */
/*     DATA IMACH(12) /       -125 / */
/*     DATA IMACH(13) /        127 / */
/*     DATA IMACH(14) /         53 / */
/*     DATA IMACH(15) /      -1021 / */
/*     DATA IMACH(16) /       1023 / */

/*     MACHINE CONSTANTS FOR THE IBM RS 6000 */

/*     DATA IMACH( 1) /          5 / */
/*     DATA IMACH( 2) /          6 / */
/*     DATA IMACH( 3) /          6 / */
/*     DATA IMACH( 4) /          0 / */
/*     DATA IMACH( 5) /         32 / */
/*     DATA IMACH( 6) /          4 / */
/*     DATA IMACH( 7) /          2 / */
/*     DATA IMACH( 8) /         31 / */
/*     DATA IMACH( 9) / 2147483647 / */
/*     DATA IMACH(10) /          2 / */
/*     DATA IMACH(11) /         24 / */
/*     DATA IMACH(12) /       -125 / */
/*     DATA IMACH(13) /        128 / */
/*     DATA IMACH(14) /         53 / */
/*     DATA IMACH(15) /      -1021 / */
/*     DATA IMACH(16) /       1024 / */

/*     MACHINE CONSTANTS FOR THE INTEL i860 */

/*     DATA IMACH( 1) /          5 / */
/*     DATA IMACH( 2) /          6 / */
/*     DATA IMACH( 3) /          6 / */
/*     DATA IMACH( 4) /          6 / */
/*     DATA IMACH( 5) /         32 / */
/*     DATA IMACH( 6) /          4 / */
/*     DATA IMACH( 7) /          2 / */
/*     DATA IMACH( 8) /         31 / */
/*     DATA IMACH( 9) / 2147483647 / */
/*     DATA IMACH(10) /          2 / */
/*     DATA IMACH(11) /         24 / */
/*     DATA IMACH(12) /       -125 / */
/*     DATA IMACH(13) /        128 / */
/*     DATA IMACH(14) /         53 / */
/*     DATA IMACH(15) /      -1021 / */
/*     DATA IMACH(16) /       1024 / */

/*     MACHINE CONSTANTS FOR THE PDP-10 (KA PROCESSOR) */

/*     DATA IMACH( 1) /          5 / */
/*     DATA IMACH( 2) /          6 / */
/*     DATA IMACH( 3) /          5 / */
/*     DATA IMACH( 4) /          6 / */
/*     DATA IMACH( 5) /         36 / */
/*     DATA IMACH( 6) /          5 / */
/*     DATA IMACH( 7) /          2 / */
/*     DATA IMACH( 8) /         35 / */
/*     DATA IMACH( 9) / "377777777777 / */
/*     DATA IMACH(10) /          2 / */
/*     DATA IMACH(11) /         27 / */
/*     DATA IMACH(12) /       -128 / */
/*     DATA IMACH(13) /        127 / */
/*     DATA IMACH(14) /         54 / */
/*     DATA IMACH(15) /       -101 / */
/*     DATA IMACH(16) /        127 / */

/*     MACHINE CONSTANTS FOR THE PDP-10 (KI PROCESSOR) */

/*     DATA IMACH( 1) /          5 / */
/*     DATA IMACH( 2) /          6 / */
/*     DATA IMACH( 3) /          5 / */
/*     DATA IMACH( 4) /          6 / */
/*     DATA IMACH( 5) /         36 / */
/*     DATA IMACH( 6) /          5 / */
/*     DATA IMACH( 7) /          2 / */
/*     DATA IMACH( 8) /         35 / */
/*     DATA IMACH( 9) / "377777777777 / */
/*     DATA IMACH(10) /          2 / */
/*     DATA IMACH(11) /         27 / */
/*     DATA IMACH(12) /       -128 / */
/*     DATA IMACH(13) /        127 / */
/*     DATA IMACH(14) /         62 / */
/*     DATA IMACH(15) /       -128 / */
/*     DATA IMACH(16) /        127 / */

/*     MACHINE CONSTANTS FOR PDP-11 FORTRAN SUPPORTING */
/*     32-BIT INTEGER ARITHMETIC. */

/*     DATA IMACH( 1) /          5 / */
/*     DATA IMACH( 2) /          6 / */
/*     DATA IMACH( 3) /          5 / */
/*     DATA IMACH( 4) /          6 / */
/*     DATA IMACH( 5) /         32 / */
/*     DATA IMACH( 6) /          4 / */
/*     DATA IMACH( 7) /          2 / */
/*     DATA IMACH( 8) /         31 / */
/*     DATA IMACH( 9) / 2147483647 / */
/*     DATA IMACH(10) /          2 / */
/*     DATA IMACH(11) /         24 / */
/*     DATA IMACH(12) /       -127 / */
/*     DATA IMACH(13) /        127 / */
/*     DATA IMACH(14) /         56 / */
/*     DATA IMACH(15) /       -127 / */
/*     DATA IMACH(16) /        127 / */

/*     MACHINE CONSTANTS FOR PDP-11 FORTRAN SUPPORTING */
/*     16-BIT INTEGER ARITHMETIC. */

/*     DATA IMACH( 1) /          5 / */
/*     DATA IMACH( 2) /          6 / */
/*     DATA IMACH( 3) /          5 / */
/*     DATA IMACH( 4) /          6 / */
/*     DATA IMACH( 5) /         16 / */
/*     DATA IMACH( 6) /          2 / */
/*     DATA IMACH( 7) /          2 / */
/*     DATA IMACH( 8) /         15 / */
/*     DATA IMACH( 9) /      32767 / */
/*     DATA IMACH(10) /          2 / */
/*     DATA IMACH(11) /         24 / */
/*     DATA IMACH(12) /       -127 / */
/*     DATA IMACH(13) /        127 / */
/*     DATA IMACH(14) /         56 / */
/*     DATA IMACH(15) /       -127 / */
/*     DATA IMACH(16) /        127 / */

/*     MACHINE CONSTANTS FOR THE SILICON GRAPHICS */

/*     DATA IMACH( 1) /          5 / */
/*     DATA IMACH( 2) /          6 / */
/*     DATA IMACH( 3) /          6 / */
/*     DATA IMACH( 4) /          6 / */
/*     DATA IMACH( 5) /         32 / */
/*     DATA IMACH( 6) /          4 / */
/*     DATA IMACH( 7) /          2 / */
/*     DATA IMACH( 8) /         31 / */
/*     DATA IMACH( 9) / 2147483647 / */
/*     DATA IMACH(10) /          2 / */
/*     DATA IMACH(11) /         24 / */
/*     DATA IMACH(12) /       -125 / */
/*     DATA IMACH(13) /        128 / */
/*     DATA IMACH(14) /         53 / */
/*     DATA IMACH(15) /      -1021 / */
/*     DATA IMACH(16) /       1024 / */

/*     MACHINE CONSTANTS FOR THE SUN */

/*     DATA IMACH( 1) /          5 / */
/*     DATA IMACH( 2) /          6 / */
/*     DATA IMACH( 3) /          6 / */
/*     DATA IMACH( 4) /          6 / */
/*     DATA IMACH( 5) /         32 / */
/*     DATA IMACH( 6) /          4 / */
/*     DATA IMACH( 7) /          2 / */
/*     DATA IMACH( 8) /         31 / */
/*     DATA IMACH( 9) / 2147483647 / */
/*     DATA IMACH(10) /          2 / */
/*     DATA IMACH(11) /         24 / */
/*     DATA IMACH(12) /       -125 / */
/*     DATA IMACH(13) /        128 / */
/*     DATA IMACH(14) /         53 / */
/*     DATA IMACH(15) /      -1021 / */
/*     DATA IMACH(16) /       1024 / */

/*     MACHINE CONSTANTS FOR THE SUN */
/*     USING THE -r8 COMPILER OPTION */

/*     DATA IMACH( 1) /          5 / */
/*     DATA IMACH( 2) /          6 / */
/*     DATA IMACH( 3) /          6 / */
/*     DATA IMACH( 4) /          6 / */
/*     DATA IMACH( 5) /         32 / */
/*     DATA IMACH( 6) /          4 / */
/*     DATA IMACH( 7) /          2 / */
/*     DATA IMACH( 8) /         31 / */
/*     DATA IMACH( 9) / 2147483647 / */
/*     DATA IMACH(10) /          2 / */
/*     DATA IMACH(11) /         53 / */
/*     DATA IMACH(12) /      -1021 / */
/*     DATA IMACH(13) /       1024 / */
/*     DATA IMACH(14) /        113 / */
/*     DATA IMACH(15) /     -16381 / */
/*     DATA IMACH(16) /      16384 / */

/*     MACHINE CONSTANTS FOR THE UNIVAC 1100 SERIES FTN COMPILER */

/*     DATA IMACH( 1) /          5 / */
/*     DATA IMACH( 2) /          6 / */
/*     DATA IMACH( 3) /          1 / */
/*     DATA IMACH( 4) /          6 / */
/*     DATA IMACH( 5) /         36 / */
/*     DATA IMACH( 6) /          4 / */
/*     DATA IMACH( 7) /          2 / */
/*     DATA IMACH( 8) /         35 / */
/*     DATA IMACH( 9) / O377777777777 / */
/*     DATA IMACH(10) /          2 / */
/*     DATA IMACH(11) /         27 / */
/*     DATA IMACH(12) /       -128 / */
/*     DATA IMACH(13) /        127 / */
/*     DATA IMACH(14) /         60 / */
/*     DATA IMACH(15) /      -1024 / */
/*     DATA IMACH(16) /       1023 / */

/*     MACHINE CONSTANTS FOR THE Z80 MICROPROCESSOR */

/*     DATA IMACH( 1) /          1 / */
/*     DATA IMACH( 2) /          1 / */
/*     DATA IMACH( 3) /          0 / */
/*     DATA IMACH( 4) /          1 / */
/*     DATA IMACH( 5) /         16 / */
/*     DATA IMACH( 6) /          2 / */
/*     DATA IMACH( 7) /          2 / */
/*     DATA IMACH( 8) /         15 / */
/*     DATA IMACH( 9) /      32767 / */
/*     DATA IMACH(10) /          2 / */
/*     DATA IMACH(11) /         24 / */
/*     DATA IMACH(12) /       -127 / */
/*     DATA IMACH(13) /        127 / */
/*     DATA IMACH(14) /         56 / */
/*     DATA IMACH(15) /       -127 / */
/*     DATA IMACH(16) /        127 / */

/* ***FIRST EXECUTABLE STATEMENT  I1MACH */
    if (*i__ < 1 || *i__ > 16) {
	goto L10;
    }

    ret_val = imach[*i__ - 1];
    return ret_val;

L10:
    //    io___3.ciunit = *output;
    //    s_wsfe(&io___3);
    //    e_wsfe();

/*     CALL FDUMP */

//    s_stop("", (ftnlen)0);
    return ret_val;
} /* i1mach_ */

#undef output
#undef imach


