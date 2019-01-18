#include<stdio.h>
#include<math.h>
#include"../Headers/vector.h"
#include"../Headers/main.h"
    
int funct_Q(double *a, double *uu)
{
   int i;
   double r;
   double n, p, u=0, v=0, w=0;
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

   a[0] = n;
   a[1] = ((K-1)*n*pow(w,2.0)+(K-1)*n*pow(v,2.0)+(K-1)*n*pow(u,2.0)+2*p)/(2*K-2);
#if dim == 1
   a[2] = n*u;
#elif dim == 2
   a[2] = n*u;
   a[3] = n*v;
#elif dim == 3 || dim == 4
   a[2] = n*u;
   a[3] = n*v;
   a[4] = n*w;
#endif

   return 0;
}
