#include<stdio.h>
#include<math.h>
#include"../Headers/matrix.h"
#include"../Headers/main.h"
    
int funct_Dm(double *a, double *uu)
{
   int i;   double r;
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

   a[0] = -(sqrt(K*n*p)-n*u)/n;
   a[1] = (n*u+sqrt(K*n*p))/n;
   a[2] = u;
     
   return 0;
}
