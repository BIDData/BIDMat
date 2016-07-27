#include "f2c.h"
#include <math.h>

__device__ integer c__1 = 1;
__device__ integer c__2 = 2;
__device__ integer c__3 = 3;
__device__ integer c__4 = 4;
__device__ integer c__5 = 5;
__device__ integer c__8 = 8;
__device__ integer c__11 = 11;
__device__ integer c__12 = 12;
__device__ integer c__13 = 13;
__device__ integer c__16 = 16;
__device__ integer c__23 = 23;

__device__ __forceinline__ double r_int(real *x) {
  return( (*x>0) ? floor(*x) : -floor(- *x) );
}

__device__ __forceinline__ double r_mod(real *x, real *y) {
  double quotient;
  if( (quotient = (double)*x / *y) >= 0)
    quotient = floor(quotient);
  else
    quotient = -floor(-quotient);
  return(*x - (*y) * quotient );
}

__device__  __forceinline__ double r_sign(real *a, real *b) {
  double x;
  x = (*a >= 0 ? *a : - *a);
  return( *b >= 0 ? x : -x);
}

__device__  __forceinline__ double pow_ri(real *ap, integer *bp) {
double pow, x;
integer n;
unsigned long u;

pow = 1;
x = *ap;
n = *bp;

if(n != 0)
        {
        if(n < 0)
                {
                n = -n;
                x = 1/x;
                }
        for(u = n; ; )
                {
                if(u & 01)
                        pow *= x;
                if(u >>= 1)
                        x *= x;
                else
                        break;
                }
        }
return(pow);
}