#include<stdio.h>
#include<math.h>
#include"../Headers/vector.h"
#include"../Headers/main.h"
    
int funct_F(double *a, double *uu)
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

   a[0] = n*u;
   a[1] = ((K-1)*n*u*pow(w,2.0)+(K-1)*n*u*pow(v,2.0)+(K-1)*n*pow(u,3.0)+2*K*p*u)/(2*K-2);
#if dim == 1
   a[2] = n*pow(u,2.0)+p;
#elif dim == 2
   a[2] = n*pow(u,2.0)+p;
   a[3] = n*u*v;
#elif dim == 3 || dim == 4
   a[2] = n*pow(u,2.0)+p;
   a[3] = n*u*v;
   a[4] = n*u*w;
#endif

   return 0;
}
