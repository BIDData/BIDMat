#include <math.h>
#include "f2c.h"

double r_int(real *x) {
  return( (*x>0) ? floor(*x) : -floor(- *x) );
}

double r_mod(real *x, real *y) {
  double quotient;
  if( (quotient = (double)*x / *y) >= 0)
    quotient = floor(quotient);
  else
    quotient = -floor(-quotient);
  return(*x - (*y) * quotient );
}

double r_sign(real *a, real *b) {
  double x;
  x = (*a >= 0 ? *a : - *a);
  return( *b >= 0 ? x : -x);
}

double pow_ri(real *ap, integer *bp) {
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
