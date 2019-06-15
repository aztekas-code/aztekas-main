#include"main.h"
    
void Prim2FluxH(double *a, double *uu)
{
   int i;
   double r;
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

   a[0] = n*w;
   a[1] = ((K-1)*n*pow(w,3.0)+((K-1)*n*pow(v,2.0)+(K-1)*n*pow(u,2.0)+2*K*p)*w)/(2*K-2);
   a[2] = n*u*w;
   a[3] = n*v*w;
   a[4] = n*pow(w,2.0)+p;
}
