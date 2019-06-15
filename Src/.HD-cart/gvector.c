#include"main.h"
    
void Prim2FluxG(double *a, double *uu)
{
   int i;
   double r;
   double n=0.0, p=0.0, u=0.0, v=0.0, w=0.0;
   n = uu[0];
   p = uu[1];

#if dim == 1
   u = uu[2];
#elif dim == 2
   u = uu[2];
   v = uu[3];
#elif dim == 3 || dim == 4
   u = uu[2];
   v = uu[3];
   w = uu[4];
#endif

   a[0] = n*v;
   a[1] = ((K-1)*n*v*pow(w,2.0)+(K-1)*n*pow(v,3.0)+((K-1)*n*pow(u,2.0)+2*K*p)*v)/(2*K-2);
#if dim == 1
   a[2] = n*u*v;
#elif dim == 2
   a[2] = n*u*v;
   a[3] = n*pow(v,2.0)+p;
#elif dim == 3 || dim == 4
   a[2] = n*u*v;
   a[3] = n*pow(v,2.0)+p;
   a[4] = n*v*w;
#endif
}
