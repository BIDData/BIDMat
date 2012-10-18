/* Copyright (c) 2012, Regents of the University of California                     */
/* All rights reserved.                                                            */

/* Redistribution and use in source and binary forms, with or without              */
/* modification, are permitted provided that the following conditions are met:     */
/*     * Redistributions of source code must retain the above copyright            */
/*       notice, this list of conditions and the following disclaimer.             */
/*     * Redistributions in binary form must reproduce the above copyright         */
/*       notice, this list of conditions and the following disclaimer in the       */
/*       documentation and/or other materials provided with the distribution.      */
/*     * Neither the name of the <organization> nor the                            */
/*       names of its contributors may be used to endorse or promote products      */
/*       derived from this software without specific prior written permission.     */

/* THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND */
/* ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED   */
/* WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE          */
/* DISCLAIMED. IN NO EVENT SHALL <COPYRIGHT HOLDER> BE LIABLE FOR ANY              */
/* DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES      */
/* (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;    */
/* LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND     */
/* ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT      */
/* (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS   */
/* SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.                    */

package edu.berkeley.bid;

public final class LAPACK {

    private LAPACK() {}

    static {
        System.loadLibrary("bidmatmkl");
    }

/*    public final static class ORDER {
        private ORDER() {}
        public final static int RowMajor=101;
        public final static int ColMajor=102;
    } */

  public static native  int sgetrf( int order, int M, int N, float [] A, int lda, int [] ipiv);
  public static native  int dgetrf( int order, int M, int N, double [] A, int lda, int [] ipiv);
  public static native  int cgetrf( int order, int M, int N, float [] A, int lda, int [] ipiv);
  public static native  int zgetrf( int order, int M, int N, double [] A, int lda, int [] ipiv);
  
  public static native  int sgetri( int order, int N, float [] A, int lda, int [] ipiv);
  public static native  int dgetri( int order, int N, double [] A, int lda, int [] ipiv);
  public static native  int cgetri( int order, int N, float [] A, int lda, int [] ipiv);
  public static native  int zgetri( int order, int N, double [] A, int lda, int [] ipiv);
  
  public static native  int sgetrs( int order, String trans, int N, int nrhs, float [] A, int lda, int [] ipiv, float [] b, int ldb);
  public static native  int dgetrs( int order, String trans, int N, int nrhs, double [] A, int lda, int [] ipiv, double [] b, int ldb);
  public static native  int cgetrs( int order, String trans, int N, int nrhs, float [] A, int lda, int [] ipiv, float [] b, int ldb);
  public static native  int zgetrs( int order, String trans, int N, int nrhs, double [] A, int lda, int [] ipiv, double [] b, int ldb);
  
  public static native  int strtrs( int order, String mdata, int n, int nrhs, float [] A, int lda, float [] b, int ldb);
  public static native  int dtrtrs( int order, String mdata, int n, int nrhs, double [] A, int lda, double [] b, int ldb);
  public static native  int ctrtrs( int order, String mdata, int n, int nrhs, float [] A, int lda, float [] b, int ldb);
  public static native  int ztrtrs( int order, String mdata, int n, int nrhs, double [] A, int lda, double [] b, int ldb);
  
  public static native  int ssteqr( int order, String compz, int n, float [] d, float [] e, float [] z, int ldz );
  public static native  int dsteqr( int order, String compz, int n, double [] d, double [] e, double [] z, int ldz );
  public static native  int csteqr( int order, String compz, int n, float [] d, float [] e, float [] z, int ldz );
  public static native  int zsteqr( int order, String compz, int n, double [] d, double [] e, double [] z, int ldz );
  
  public static native  int ssytrd( int order, String uplo, int n, float [] a, int lda, float [] d, float [] e, float [] tau );
  public static native  int dsytrd( int order, String uplo, int n, double [] a, int lda, double [] d, double [] e, double [] tau );
  
  public static native  int sorgtr( int order, String uplo, int n, float [] a, int lda, float [] tau );
  public static native  int dorgtr( int order, String uplo, int n, double [] a, int lda, double [] tau );
  
  public static native  int sstedc( int order, String compz, int n, float [] d, float [] e, float [] z, int ldz );
  public static native  int dstedc( int order, String compz, int n, double [] d, double [] e, double [] z, int ldz );
  
  public static native  int ssyevd( int order, String jobz, String uplo, int n, float [] a, int lda, float [] w );
  public static native  int dsyevd( int order, String jobz, String uplo, int n, double [] a, int lda, double [] w );
  
  public static native  int spotrf( int order, String uplo, int n, float [] a, int lda);
  public static native  int dpotrf( int order, String uplo, int n, double [] a, int lda);
  public static native  int cpotrf( int order, String uplo, int n, float [] a, int lda); 
  public static native  int zpotrf( int order, String uplo, int n, double [] a, int lda); 

  public static native int sgebal(int matrix_order, String job, int n, float [] a, int lda, int [] ilo, int [] ihi, float [] scale);
  public static native int dgebal(int matrix_order, String job, int n, double [] a, int lda, int [] ilo, int [] ihi, double [] scale);
  public static native int cgebal(int matrix_order, String job, int n, float [] a, int lda, int [] ilo, int [] ihi, float [] scale);
  public static native int zgebal(int matrix_order, String job, int n, double [] a, int lda, int [] ilo, int [] ihi, double [] scale);
  
  public static native int cunghr(int matrix_order, int n, int ilo, int ihi, float [] a, int lda, float [] tau);
  public static native int zunghr(int matrix_order, int n, int ilo, int ihi, double [] a, int lda, double [] tau);

  public static native int strevc(int matrix_order, String side, String howmny, int [] select, int n, float [] t, int ldt, float [] vl, int ldvl, float [] vr, int ldvr, int mm, int [] m);
  public static native int dtrevc(int matrix_order, String side, String howmny, int [] select, int n, double [] t, int ldt, double [] vl, int ldvl, double [] vr, int ldvr, int mm, int [] m);
  public static native int ctrevc(int matrix_order, String side, String howmny, int [] select, int n, float [] t, int ldt, float [] vl, int ldvl, float [] vr, int ldvr, int mm, int [] m);
  public static native int ztrevc(int matrix_order, String side, String howmny, int [] select, int n, double [] t, int ldt, double [] vl, int ldvl, double [] vr, int ldvr, int mm, int [] m);

  public static native int sgehrd(int matrix_order, int n, int ilo, int ihi, float [] a, int lda, float [] tau);
  public static native int dgehrd(int matrix_order, int n, int ilo, int ihi, double [] a, int lda, double [] tau);
  public static native int cgehrd(int matrix_order, int n, int ilo, int ihi, float [] a, int lda, float [] tau);
  public static native int zgehrd(int matrix_order, int n, int ilo, int ihi, double [] a, int lda, double [] tau);

  public static native int shseqr(int matrix_order, String job, String compz, int n, int ilo, int ihi, float [] h, int ldh, float [] wr, float [] wi, float [] z, int ldz);
  public static native int dhseqr(int matrix_order, String job, String compz, int n, int ilo, int ihi, double [] h, int ldh, double [] wr, double [] wi, double [] z, int ldz);
  public static native int chseqr(int matrix_order, String job, String compz, int n, int ilo, int ihi, float [] h, int ldh, float [] w, float [] z, int ldz);
  public static native int zhseqr(int matrix_order, String job, String compz, int n, int ilo, int ihi, double [] h, int ldh, double [] w, double [] z, int ldz);

  public static native int sgebak(int matrix_order, String job, String side, int n, int ilo, int ihi, float [] scale, int m, float [] v, int ldv);
  public static native int dgebak(int matrix_order, String job, String side, int n, int ilo, int ihi, double [] scale, int m, double [] v, int ldv);
  public static native int cgebak(int matrix_order, String job, String side, int n, int ilo, int ihi, float [] scale, int m, float [] v, int ldv);
  public static native int zgebak(int matrix_order, String job, String side, int n, int ilo, int ihi, double [] scale, int m, double [] v, int ldv);

  public static native int sgeqrf(int matrix_order, int m, int n, float [] a, int lda, float [] tau);
  public static native int dgeqrf(int matrix_order, int m, int n, double [] a, int lda, double [] tau);
  public static native int cgeqrf(int matrix_order, int m, int n, float [] a, int lda, float [] tau);
  public static native int zgeqrf(int matrix_order, int m, int n, double [] a, int lda, double [] tau);

  public static native int sgeqp3(int matrix_order, int m, int n, float [] a, int lda, int [] jpvt, float [] tau);
  public static native int dgeqp3(int matrix_order, int m, int n, double [] a, int lda, int [] jpvt, double [] tau);
  public static native int cgeqp3(int matrix_order, int m, int n, float [] a, int lda, int [] jpvt, float [] tau);
  public static native int zgeqp3(int matrix_order, int m, int n, double [] a, int lda, int [] jpvt, double [] tau);

  public static native int sorgqr(int matrix_order, int m, int n, int k, float [] a, int lda, float [] tau);
  public static native int dorgqr(int matrix_order, int m, int n, int k, double [] a, int lda, double [] tau);

  public static native int cungqr(int matrix_order, int m, int n, int k, float [] a, int lda, float [] tau);
  public static native int zungqr(int matrix_order, int m, int n, int k, double [] a, int lda, double [] tau);
}