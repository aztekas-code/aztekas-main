#include<stdio.h>
#include<math.h>
#include"../Headers/matrix.h"
#include"../Headers/main.h"
    
int funct_Dn(double *a, double *uu)
{
   int i;   double r;
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

   a[0] = -(sqrt(K*n*p)-n*v)/n;
   a[1] = (n*v+sqrt(K*n*p))/n;
   a[2] = v;

     
   return 0;
}
