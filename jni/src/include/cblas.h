// Android CBLAS header
// Points to defined CBLAS implementation and back-fills missing/different
// fuctions in the implementations that require it.

#ifndef BIDMAT_CBLAS_WRAPPER_H
#define BIDMAT_CBLAS_WRAPPER_H

#ifdef QSML

// Ensure parity when compiling as C or C++
#ifndef __cplusplus
#include <complex.h>

// weird bug, won't compile C without this
#define complex __complex__

#include <qblas_cblas_common.h>
typedef enum CBLAS_ORDER CBLAS_ORDER;
typedef enum CBLAS_TRANSPOSE CBLAS_TRANSPOSE;
typedef enum CBLAS_UPLO CBLAS_UPLO;
typedef enum CBLAS_DIAG CBLAS_DIAG;
typedef enum CBLAS_SIDE CBLAS_SIDE;
#endif

#include <qblas_cblas.h>

#elif OpenBLAS
#include <openblas_cblas.h>
#endif

#endif // BIDMAT_CBLAS_WRAPPER_H
