#include<stdio.h>
#include<math.h>
#include"../Headers/vector.h"
#include"../Headers/main.h"
    
int funct_S(double *a, double *uu)
{
   int i;
   double r;
   double n, p, u=0, v=0, w=0;
   n = uu[0];
   p = uu[1];
   u = uu[2];
   r = sqrt(x1*x1 + x2*x2);
   if(dim >= 2){v = uu[3];}
   if(dim == 3){w = uu[4];}

   a[0] = 0;
   a[1] = 0;
   a[2] = 0;
   a[3] = 0;
   a[4] = 0;

   return 0;
}
