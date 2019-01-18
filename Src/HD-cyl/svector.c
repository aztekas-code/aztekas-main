#include<stdio.h>
#include<math.h>
#include"../Headers/vector.h"
#include"../Headers/main.h"
    
int funct_S(double *a, double *uu)
{
   int i;
   double R = x1;
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

   a[0] = 0;
   a[1] = 0;
   a[2] = (n*pow(w,2.0)+p)/x1;
   a[3] = 0;
   a[4] = -n*u*w/x1;

   return 0;
}
