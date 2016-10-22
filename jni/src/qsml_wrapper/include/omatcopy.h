#ifndef OMATCOPY_H
#define OMATCOPY_H

#include <qsml_types.h>

int cblas_domatcopy(CBLAS_ORDER corder, CBLAS_TRANSPOSE ctrans, qsml_long rows, qsml_long cols, double alpha, double *a, qsml_long lda, double *b, qsml_long ldb);
int cblas_somatcopy(CBLAS_ORDER corder, CBLAS_TRANSPOSE ctrans, qsml_long rows, qsml_long cols, float alpha, float *a, qsml_long lda, float *b, qsml_long ldb);

int somatcopy_rn(qsml_long rows, qsml_long cols, float alpha, float *a, qsml_long lda, float *b, qsml_long ldb);
int somatcopy_rt(qsml_long rows, qsml_long cols, float alpha, float *a, qsml_long lda, float *b, qsml_long ldb);
int somatcopy_ct(qsml_long rows, qsml_long cols, float alpha, float *a, qsml_long lda, float *b, qsml_long ldb);
int somatcopy_cn(qsml_long rows, qsml_long cols, float alpha, float *a, qsml_long lda, float *b, qsml_long ldb);


int domatcopy_rn(qsml_long rows, qsml_long cols, double alpha, double *a, qsml_long lda, double *b, qsml_long ldb);
int domatcopy_rt(qsml_long rows, qsml_long cols, double alpha, double *a, qsml_long lda, double *b, qsml_long ldb);
int domatcopy_ct(qsml_long rows, qsml_long cols, double alpha, double *a, qsml_long lda, double *b, qsml_long ldb);
int domatcopy_cn(qsml_long rows, qsml_long cols, double alpha, double *a, qsml_long lda, double *b, qsml_long ldb);


#endif 
