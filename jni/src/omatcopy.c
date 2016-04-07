/***************************************************************************
Copyright (c) 2013, The OpenBLAS Project
All rights reserved.
Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are
met:
1. Redistributions of source code must retain the above copyright
notice, this list of conditions and the following disclaimer.
2. Redistributions in binary form must reproduce the above copyright
notice, this list of conditions and the following disclaimer in
the documentation and/or other materials provided with the
distribution.
3. Neither the name of the OpenBLAS project nor the names of
its contributors may be used to endorse or promote products
derived from this software without specific prior written permission.
THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE OPENBLAS PROJECT OR CONTRIBUTORS BE
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE
USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*****************************************************************************/

// #include <omatcopy.h>

#include <cblas.h>

// #include "common.h"

/*****************************************************
 * 2014/06/09 Saar
 *
 * Order rowMajor
 * No Trans
 *
******************************************************/


int cblas_domatcopy (CBLAS_ORDER corder, CBLAS_TRANSPOSE ctrans, qsml_long rows, qsml_long cols, double alpha, double *a, qsml_long lda, double *b, qsml_long ldb){
	if (corder == CblasRowMajor && ctrans == CblasTrans) {
		return domatcopy_rt(rows, cols, alpha, a, lda, b, ldb);
	} else if (corder == CblasColMajor && ctrans == CblasTrans) {
		return domatcopy_ct(rows, cols, alpha, a, lda, b, ldb);
	} else if (corder == CblasRowMajor && ctrans == CblasNoTrans) {
		return domatcopy_rn(rows, cols, alpha, a, lda, b, ldb);
	} else if (corder == CblasColMajor && ctrans == CblasNoTrans) {
		return domatcopy_cn(rows, cols, alpha, a, lda, b, ldb);
	}
}

int cblas_somatcopy (CBLAS_ORDER corder, CBLAS_TRANSPOSE ctrans, qsml_long rows, qsml_long cols, float alpha, float *a, qsml_long lda, float *b, qsml_long ldb){
	if (corder == CblasRowMajor && ctrans == CblasTrans) {
		return somatcopy_rt(rows, cols, alpha, a, lda, b, ldb);
	} else if (corder == CblasColMajor && ctrans == CblasTrans) {
		return somatcopy_ct(rows, cols, alpha, a, lda, b, ldb);
	} else if (corder == CblasRowMajor && ctrans == CblasNoTrans) {
		return somatcopy_rn(rows, cols, alpha, a, lda, b, ldb);
	} else if (corder == CblasColMajor && ctrans == CblasNoTrans) {
		return somatcopy_cn(rows, cols, alpha, a, lda, b, ldb);
	}
}

int somatcopy_rn(qsml_long rows, qsml_long cols, float alpha, float *a, qsml_long lda, float *b, qsml_long ldb)
{
	qsml_long i,j;
	float *aptr,*bptr;

	if ( rows <= 0     )  return(0);
	if ( cols <= 0     )  return(0);

	aptr = a;
	bptr = b;

	if ( alpha == 0.0 )
	{
		for ( i=0; i<rows ; i++ )
		{
			for(j=0; j<cols; j++)
			{
				bptr[j] = 0.0;
			}
			bptr += ldb;
		}
		return(0);
	}

	if ( alpha == 1.0 )
	{
		for ( i=0; i<rows ; i++ )
		{
			for(j=0; j<cols; j++)
			{
				bptr[j] = aptr[j];
			}
			aptr += lda;
			bptr += ldb;
		}
		return(0);
	}

	for ( i=0; i<rows ; i++ )
	{
		for(j=0; j<cols; j++)
		{
			bptr[j] = alpha * aptr[j];
		}
		aptr += lda;
		bptr += ldb;
	}

	return(0);

}

int somatcopy_rt(qsml_long rows, qsml_long cols, float alpha, float *a, qsml_long lda, float *b, qsml_long ldb)
{
	qsml_long i,j;
	float *aptr,*bptr;

	if ( rows <= 0     )  return(0);
	if ( cols <= 0     )  return(0);

	aptr = a;

	for ( i=0; i<rows ; i++ )
	{
		bptr = &b[i];
		for(j=0; j<cols; j++)
		{
			bptr[j*ldb] = alpha * aptr[j];
		}
		aptr += lda;
	}

	return(0);

}

int somatcopy_ct(qsml_long rows, qsml_long cols, float alpha, float *a, qsml_long lda, float *b, qsml_long ldb)
{
	qsml_long i,j;
	float *aptr,*bptr;

	if ( rows <= 0     )  return(0);
	if ( cols <= 0     )  return(0);

	aptr = a;

	if ( alpha == 0.0 )
	{
		for ( i=0; i<cols ; i++ )
		{
			bptr = &b[i];
			for(j=0; j<rows; j++)
			{
				bptr[j*ldb] = 0.0;
			}
		}
		return(0);
	}

	if ( alpha == 1.0 )
	{
		for ( i=0; i<cols ; i++ )
		{
			bptr = &b[i];
			for(j=0; j<rows; j++)
			{
				bptr[j*ldb] = aptr[j];
			}
			aptr += lda;
		}
		return(0);
	}

	for ( i=0; i<cols ; i++ )
	{
		bptr = &b[i];
		for(j=0; j<rows; j++)
		{
			bptr[j*ldb] = alpha * aptr[j];
		}
		aptr += lda;
	}

	return(0);

}

int somatcopy_cn(qsml_long rows, qsml_long cols, float alpha, float *a, qsml_long lda, float *b, qsml_long ldb)
{
	qsml_long i,j;
	float *aptr,*bptr;

	if ( rows <= 0     )  return(0);
	if ( cols <= 0     )  return(0);

	aptr = a;
	bptr = b;

	if ( alpha == 0.0 )
	{
		for ( i=0; i<cols ; i++ )
		{
			for(j=0; j<rows; j++)
			{
				bptr[j] = 0.0;
			}
			bptr += ldb;
		}
		return(0);
	}

	if ( alpha == 1.0 )
	{
		for ( i=0; i<cols ; i++ )
		{
			for(j=0; j<rows; j++)
			{
				bptr[j] = aptr[j];
			}
			aptr += lda;
			bptr += ldb;
		}
		return(0);
	}

	for ( i=0; i<cols ; i++ )
	{
		for(j=0; j<rows; j++)
		{
			bptr[j] = alpha * aptr[j];
		}
		aptr += lda;
		bptr += ldb;
	}

	return(0);
}




















int domatcopy_rn(qsml_long rows, qsml_long cols, double alpha, double *a, qsml_long lda, double *b, qsml_long ldb)
{
	qsml_long i,j;
	double *aptr,*bptr;

	if ( rows <= 0     )  return(0);
	if ( cols <= 0     )  return(0);

	aptr = a;
	bptr = b;

	if ( alpha == 0.0 )
	{
		for ( i=0; i<rows ; i++ )
		{
			for(j=0; j<cols; j++)
			{
				bptr[j] = 0.0;
			}
			bptr += ldb;
		}
		return(0);
	}

	if ( alpha == 1.0 )
	{
		for ( i=0; i<rows ; i++ )
		{
			for(j=0; j<cols; j++)
			{
				bptr[j] = aptr[j];
			}
			aptr += lda;
			bptr += ldb;
		}
		return(0);
	}

	for ( i=0; i<rows ; i++ )
	{
		for(j=0; j<cols; j++)
		{
			bptr[j] = alpha * aptr[j];
		}
		aptr += lda;
		bptr += ldb;
	}

	return(0);

}

int domatcopy_rt(qsml_long rows, qsml_long cols, double alpha, double *a, qsml_long lda, double *b, qsml_long ldb)
{
	qsml_long i,j;
	double *aptr,*bptr;

	if ( rows <= 0     )  return(0);
	if ( cols <= 0     )  return(0);

	aptr = a;

	for ( i=0; i<rows ; i++ )
	{
		bptr = &b[i];
		for(j=0; j<cols; j++)
		{
			bptr[j*ldb] = alpha * aptr[j];
		}
		aptr += lda;
	}

	return(0);

}

int domatcopy_ct(qsml_long rows, qsml_long cols, double alpha, double *a, qsml_long lda, double *b, qsml_long ldb)
{
	qsml_long i,j;
	double *aptr,*bptr;

	if ( rows <= 0     )  return(0);
	if ( cols <= 0     )  return(0);

	aptr = a;

	if ( alpha == 0.0 )
	{
		for ( i=0; i<cols ; i++ )
		{
			bptr = &b[i];
			for(j=0; j<rows; j++)
			{
				bptr[j*ldb] = 0.0;
			}
		}
		return(0);
	}

	if ( alpha == 1.0 )
	{
		for ( i=0; i<cols ; i++ )
		{
			bptr = &b[i];
			for(j=0; j<rows; j++)
			{
				bptr[j*ldb] = aptr[j];
			}
			aptr += lda;
		}
		return(0);
	}

	for ( i=0; i<cols ; i++ )
	{
		bptr = &b[i];
		for(j=0; j<rows; j++)
		{
			bptr[j*ldb] = alpha * aptr[j];
		}
		aptr += lda;
	}

	return(0);

}

int domatcopy_cn(qsml_long rows, qsml_long cols, double alpha, double *a, qsml_long lda, double *b, qsml_long ldb)
{
	qsml_long i,j;
	double *aptr,*bptr;

	if ( rows <= 0     )  return(0);
	if ( cols <= 0     )  return(0);

	aptr = a;
	bptr = b;

	if ( alpha == 0.0 )
	{
		for ( i=0; i<cols ; i++ )
		{
			for(j=0; j<rows; j++)
			{
				bptr[j] = 0.0;
			}
			bptr += ldb;
		}
		return(0);
	}

	if ( alpha == 1.0 )
	{
		for ( i=0; i<cols ; i++ )
		{
			for(j=0; j<rows; j++)
			{
				bptr[j] = aptr[j];
			}
			aptr += lda;
			bptr += ldb;
		}
		return(0);
	}

	for ( i=0; i<cols ; i++ )
	{
		for(j=0; j<rows; j++)
		{
			bptr[j] = alpha * aptr[j];
		}
		aptr += lda;
		bptr += ldb;
	}

	return(0);

}