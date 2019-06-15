#include"main.h"
    
void Source_Terms(double *a, double *uu)
{
   int i;
   double n, p, u=0, v=0, w=0;
   n = uu[0];
   p = uu[1];

#if DIM == 1
   u = uu[2];
#elif DIM == 2
   u = uu[2];
   v = uu[3];
#elif DIM == 3 || DIM == 4
   u = uu[2];
   v = uu[3];
   w = uu[4];
#endif

   double E = 0.5*n*(u*u + v*v + w*w) + p/(K-1);

   a[0] = 0.0;
   a[1] = 0.0;
   a[2] = 0.0;
   a[3] = 0.0;
   a[4] = 0.0;
}
